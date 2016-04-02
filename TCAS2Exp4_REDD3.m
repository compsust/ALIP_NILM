clear all
close all
numSampVA = 376150;     % number of samples
Ni = 1;
rngS = [Ni,1,Ni-1+numSampVA,14]; % range
%rngS = [1,1,numSampVA,12]; % range
datafile = ['REDD\REDDhouse3_3sec_VA.csv'; ];

%% generating energy meter output y
numApl = 7; % number of appliances
data = csvread(datafile,Ni,1,rngS);
dataDsamp = zeros(14,numSampVA);
for i = 1:14;
    dataDsamp(i,:) = data(:,i);
end
dataVA = zeros(numApl,numSampVA/10);
s=[5 7 9 10 12 8 11];
for i = 1:numApl;
    %figure
    dataVA(i,1:numSampVA/10) = downsample(dataDsamp(s(i),:),10); 
    %plot(1:numSampVA/10,dataVA(i,1:numSampVA/10));
    %axis('tight')
    %pause
    %close all
end
numSampVA=numSampVA/10;
%break

numState=20;                
y = sum(dataVA(:,:),1);
numTrial=1;
dataEsVa = zeros(numApl,numSampVA,numTrial);
acc = zeros(2,numTrial);
ratingVA = [1; 83; 115; 408;
    1;210; 737;
    295;4750;
    2;165;1760;
    12;945;1250;1600;
    5;700;
    1;8];
ac(1,:)= zeros(1,numState); ac(1,4)=1; ac(1,6)=1; ac(1,8)=1;ac(1,11)=1;
ac(1,14)=1;%
ac(2,:)= zeros(1,numState); ac(1,5)=1; ac(1,7)=1; ac(1,8)=1;ac(1,11)=1;
ac(3,:)=zeros(1,numState); ac(3,6)=1; ac(3,16)=1;
ac(4,:)=zeros(1,numState); ac(4,7)=1; ac(4,18)=1;
ac(5,:)=zeros(1,numState); ac(5,13)=1; ac(5,20)=1;

h = ones(numApl,1);

F = zeros(numApl,numState);
F(1,1:4)=ones(1,4);
F(2,5:7)=ones(1,3);
F(3,8:9)=ones(1,2);
F(4,10:12)=ones(1,3);
F(5,13:16)=ones(1,4);
F(6,17:18)=ones(1,2);
F(7,19:20)=ones(1,2);
u1 = ones(4,1);

%% Integer LP initialization begins
c = [1; zeros(numState,1)];
b = [1; zeros(numState,1)];
Fh=F'*h;
a = [0; (ratingVA.*Fh)];
B1 = [0 F(3,:); 0 F(5,:);0 F(7,:);];
B2 = [0 F(3,:); 0 F(5,:);0 F(7,:); 0 ac(1,:); 0 ac(2,:);
    0 ac(3,:); 0 ac(4,:); 0 ac(5,:)];
Beq = [zeros(4,1) [F(1,:); F(2,:); F(4,:); F(6,:)]];
A2 = [-(a+b)'; (a-b)'; B2];
A1 = [-(a+b)'; (a-b)'; B1; Beq];
intcon = 2:numState+1;
lb = zeros(numState+1,1); ub = [inf; ones(numState,1)];
cpu_IP = 0; cpu_ALIP=0;
for trial = 1:numTrial;
    f = ones(numState,numSampVA,2);
    zf = zeros(numApl,numSampVA,2);
    for i = 2:numSampVA;
        options=optimoptions(@intlinprog,'display','off');
        tm = cputime;
        eb = [-y(1,i);  y(1,i); 1;1;1;1;1;1;1];
        zhat = intlinprog(c,intcon,A1,eb,[],[],lb,ub,options);
        f(:,i,1)  = round(zhat(2:end,1));%
        zf(:,i,1) = F*(f(:,i,1).*ratingVA);
        elt = cputime-tm;
        cpu_IP = cpu_IP+elt;
        
        tm = cputime;
        eb = [-y(1,i);  y(1,i); 1;1;1;1;1;1;1;1];
        zhat = intlinprog(c,intcon,A2,eb,Beq,u1,lb,ub,options);
        f(:,i,2)  = round(zhat(2:end,1));%
        zf(:,i,2) = F*(f(:,i,2).*ratingVA);
        if i>3 && zf(2,i-3,2)==210 && median(zf(2,i-3:i,2))==1;
            zf(2,i-3,2) = 1;
        end
        if i>15&& zf(3,i,2)==0 && median(zf(3,i-15:i,2))==4750;
            zf(3,i,2) = 295;
        end
        if i>5 && zf(5,i-5,2)==1250 && median(zf(5,i-5:i,2))==0;
            zf(5,i-5,2) = 0;
        end
        if i>15 && zf(5,i-15,2)==12 && median(zf(5,i-15:i,2))==0;
            zf(5,i-15,2) = 0;
        end
        if zf(7,i-1,2)==0 && zf(7,i,2)==8;
            zf(7,i,2) = 1;
        end
        if i>35 && zf(7,i-35,2)==8 && median(zf(7,i-35:i,2))==1;
            zf(7,i-35,2) = 1;
        end
        elt = cputime-tm;
        cpu_ALIP = cpu_ALIP+elt;
    end
    for i=2:numSampVA
        tm = cputime;
        zft=[];
        zfn=[]; zfx=[];
        if zf(1,i,2)==115;
            zft = 115;
            zfn = 115; 
            zfx = 130;
        end
        if zf(2,i,2)==737;
            zft = [zft 737];
            zfn = [zfn 737];
            zfx = [zfx 744];            
        end
        if zf(5,i,2)==1250;
            zft = [zft 1250];
            zfn = [zfn 1250];
            zfx = [zfx 1300];            
        end
        if zf(5,i,2)==1600;
            zft = [zft 1600];
            zfn = [zfn 1580];
            zfx = [zfx 1600];            
        end
        if zf(6,i,2)==700;
            zft = [zft 700];
            zfn = [zfn 700];
            zfx = [zfx 712];            
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
            for k = 1:6;
                if zft(j)==zf(k,i,2);
                    zf(k,i,2) = zhat(j+1,1);
                end
            end
        end
        elt = cputime-tm;
        cpu_ALIP = cpu_ALIP+elt;
    end
    for j = 1:2;
        dataEsVa(:,:,trial,j) = zf(:,:,j);
        e =dataVA(:,:,trial)-zf(:,:,j);
        acc(j,trial)=1-sum(sum(abs(e),2),1)/(2*sum(abs(y(trial,:))));
    end
end
cpu_IP/numSampVA
cpu_ALIP/numSampVA
IP_ACC = acc(1,1)
ALIP_ACC = acc(2,1)

e1 =dataVA(:,:,1:numTrial)-dataEsVa(:,:,1,1);
e2 =dataVA(:,:,1:numTrial)-dataEsVa(:,:,1,2);
 
 for trial=1:numTrial;   
    for i=1:7;    
        Sacc(i,trial,1)=1-sum(abs(e1(i,:,trial)),1)/(2*sum(abs(dataVA(i,:,trial)),1));
        Sacc(i,trial,2)=1-sum(abs(e2(i,:,trial)),1)/(2*sum(abs(dataVA(i,:,trial)),1));  
    end
 end
 IP_AC = Sacc(:,1,1)
 ALIP_AC = Sacc(:,1,2)