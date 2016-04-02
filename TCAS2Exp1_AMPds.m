load 'ElectLett_exp1.mat';
dataVA = dataxExp1;
y = datayExp1;
[numApl, numSampVA,numTrial] =size(dataVA);
%numTrial = 2;
dataEsVa = zeros(numApl,numSampVA,numTrial);
acc = zeros(2,numTrial);
numState = 13;                %total states of the 3 appliances
ratingVA = [4620; 502; 29 ;
    10; 130; 445;1019;
    51; 24; 1823; 2400;
    9;22;];
minVA = [4550; 502; 
    130;
     1800; 2400;
    ];
maxVA = [4660; 520;  
    150; 
    1900; 2450;
    ];

F = zeros(numApl,numState);
F(1,1:3)=ones(1,3);
F(2,4:7)=ones(1,4);
F(3,8:11)=ones(1,4);
F(4,12:13)=ones(1,2);
up = ones(4,1);
h = ones(numApl,1);

%% Integer LP initialization begins
c = [1; zeros(numState,1)];
b = [1; zeros(numState,1)];
Fh=F'*h;
a = [0; (ratingVA.*Fh)];
%B = [zeros(2,1) ac];% [0 0 1 0 1 0 1 zeros(1,m-7) 1]];
Beq = [zeros(numApl,1) [F(1,:); F(2,:); F(3,:); F(4,:)]];
Ap = [-(a+b)'; (a-b)';];% B];
Ac = [-(a+b)'; (a-b)'; Beq];
intcon = 2:numState+1;
lb = zeros(numState+1,1); ub = [inf; ones(numState,1)];

for trial = 1:numTrial;
    trial
    f = zeros(numState,numSampVA,2);
    zf = zeros(numApl,numSampVA,2);
    for i = 2:numSampVA;
        options=optimoptions(@intlinprog,'display','off');
        ep = [-y(trial,i);  y(trial,i);];
        ec = [ep; up];
        
        zhat = intlinprog(c,intcon,Ac,ec,[],[],lb,ub,options);
        f(:,i,1)  = round(zhat(2:end,1));%
        zf(:,i,1) = F*(f(:,i,1).*ratingVA); 
        
        zhat = intlinprog(c,intcon,Ap,ep,Beq,up,lb,ub,options);
        f(:,i,2)  = round(zhat(2:end,1));%
        zf(:,i,2) = F*(f(:,i,2).*ratingVA);
        
        if zf(1,i-1,2)==29 && zf(1,i,2)==502;
            zf(1,i,2)= 29;
        end
        if zf(2,i-1,2)==10 && zf(2,i,2)==445;
            zf(2,i,2) = 130;
        end
        if zf(2,i-1,2)==10 && zf(2,i,2)==1019;
            zf(2,i,2) = 10;
        end
        if zf(2,i-1,2)==130 && zf(2,i,2)==1019;
            zf(2,i,2) = 130;
        end
        if zf(3,i-1,2) ==51 && zf(3,i,2)==24;
            zf(3,i,2) = 51;
        end
        if zf(3,i-1,2)==1823 && zf(3,i,2) == 51;
            zf(3,i,2) = 24;
        end
        if zf(3,i-1,2)==2400 && median(zf(3,i-1:i,2))<2400;
            zf(3,i-1:i,2) = 1823;
        end
        if i>45 && median(zf(4,i-44:i,2))<22;
            zf(4,i-44,2)=9;
        end
    end
    for i = 2:numSampVA;
        p2=[]; zft=[];
        zfn=[]; zfx=[];
        if zf(1,i,2)==4620;
            zft = 4620;
            zfn = 4550; 
            zfx = 4660;
            p2(1)=1;
        end
        if zf(1,i,2)==502;
            zft = [zft 502];
            zfn = [zfn 502];
            zfx = [zfx 520];
            p2(2)=2;
        end
        if zf(2,i,2)==130;
            zft = [zft 130];
            zfn = [zfn 127];
            zfx = [zfx 150];
            p2(3)=5;
        end
        if zf(3,i,2)==1823;
            zft = [zft 1823];
            zfn = [zfn 1800];
            zfx = [zfx 1900];
            p2(4)=10;
        end
        if zf(3,i,2)==2400;
            zft = [zft 2400];
            zfn = [zfn 2400];
            zfx = [zfx 2450];
            p2(5)=11;
        end
        hl = ones(1,length(zft));
        tf = isempty(zft);
        if tf==0;
            cy=y(trial,i)-(h'*zf(:,i,2)-hl*zft');
            htilde = [0; ones(length(zft),1)];
            uh = [1; zeros(length(zft),1)];
            A = [-(htilde+uh)'; (htilde-uh)';
                zeros(length(zft),1) -eye(length(zft));
                zeros(length(zft),1) eye(length(zft))];
            cr = [1 zeros(1,length(zft))];
            options = optimoptions(@linprog,'display','off');
            ep = [-cy;  cy;-zfn'; zfx'];
            zhat = linprog(cr,A,ep,[],[],[],[],[],options);
        end
        for j=1:length(zft);
            for k = 1:3;
                if zft(j)==zf(k,i,2);
                    zf(k,i,2) = zhat(j+1,1);
                end
            end
        end
    end
    for j = 1:2;
        dataEsVa(:,:,trial,j) = zf(:,:,j);
        e =dataVA(:,:,trial)-zf(:,:,j);
        acc(j,trial)=1-sum(sum(abs(e),2),1)/(2*sum(abs(y(trial,:))));
    end    
end
e1 =dataVA(:,:,1:numTrial)-dataEsVa(:,:,1:numTrial,1);
e2 =dataVA(:,:,1:numTrial)-dataEsVa(:,:,1:numTrial,2);
 for trial=1:numTrial;   
    for i=1:4;    
        Sacc(i,trial,1)=1-sum(abs(e1(i,:,trial)),1)/(2*sum(abs(dataVA(i,:,trial)),1));
        Sacc(i,trial,2)=1-sum(abs(e2(i,:,trial)),1)/(2*sum(abs(dataVA(i,:,trial)),1));
    end
 end
figure;
subplot(2,2,1);
plot(1:numTrial,Sacc(1,1:numTrial,1),'-r','linewidth',2,'MarkerSize',10);
hold on
plot(1:numTrial,Sacc(1,1:numTrial,2),'-b','linewidth',2,'MarkerSize',10);
%title('CDE')
%legend('IP', 'ALIP')
xlabel('Number of blocks');
ylabel('AC, CDE');
axis('tight')
grid on
subplot(2,2,2);
plot(1:numTrial,Sacc(2,1:numTrial,1),'-r','linewidth',2,'MarkerSize',10);
hold on
plot(1:numTrial,Sacc(2,1:numTrial,2),'-b','linewidth',2,'MarkerSize',10);
%title( 'FRG')
%legend('IP', 'ALIP')
xlabel('Number of blocks');
ylabel('AC, FRG');
axis('tight')
grid on
subplot(2,2,3);
plot(1:numTrial,Sacc(4,1:numTrial,1),'-r','linewidth',2,'MarkerSize',10);
hold on
plot(1:numTrial,Sacc(4,1:numTrial,2),'-b','linewidth',2,'MarkerSize',10)
%title('B1E')
%legend('IP', 'ALIP')
xlabel('Number of blocks');
ylabel('AC, B1E');
axis('tight')
grid on

subplot(2,2,4);
plot(1:numTrial,acc(1,1:numTrial),'-r','linewidth',2,'MarkerSize',5);
hold on
plot(1:numTrial,acc(2,1:numTrial),'-b','linewidth',2,'MarkerSize',5);
legend('IP', 'ALIP')
axis('tight')
xlabel('Number of blocks');
ylabel('ACC');
grid on

figure
plot(1:numTrial,Sacc(3,1:numTrial,1),'-r','linewidth',1.5,'MarkerSize',10);
hold on
plot(1:numTrial,Sacc(3,1:numTrial,2),'-b','linewidth',1.5,'MarkerSize',10);
title('HPE')
legend('IP', 'PIP')
xlabel('Number of blocks');
ylabel('AC');
axis('tight')
grid on

shdataVA = reshape(dataVA(:,:,1:numTrial),[4,5040*numTrial]);
shdataEsVA(:,:,1) = reshape(dataEsVa(:,:,1:numTrial,1),[4,5040*numTrial]);
shdataEsVA(:,:,2) = reshape(dataEsVa(:,:,1:numTrial,2),[4,5040*numTrial]);
e1 =shdataVA-shdataEsVA(:,:,1);
e2 =shdataVA-shdataEsVA(:,:,2);
for i=1:4;    
        Sacc(i,1)=1-sum(abs(e1(i,:)),1)/(2*sum(abs(shdataVA(i,:)),1));
        Sacc(i,2)=1-sum(abs(e2(i,:)),1)/(2*sum(abs(shdataVA(i,:)),1));
end

Nacc(1,1)=1-sum(sum(abs(e1),2),1)/(2*sum(sum(abs(shdataVA),2),1));
Nacc(2,1)=1-sum(sum(abs(e2),2),1)/(2*sum(sum(abs(shdataVA),2),1));
IP_AC = Sacc(:,1) 
ALIP_AC=Sacc(:,2)
IP_ACC = Nacc(1,1)
ALIP_ACC = Nacc(2,1)