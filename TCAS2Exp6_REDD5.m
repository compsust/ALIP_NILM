close all
clear all
numSampVA = 77450;     % number of samples
Ni = 1;
rngS = [Ni,1,Ni-1+numSampVA,16]; % range
%rngS = [1,1,numSampVA,12]; % range
datafile = ['REDD\REDDhouse5_3sec_VA.csv'; ];

%% generating energy meter output y
numApl = 6; % number of appliances
data = csvread(datafile,Ni,1,rngS);
dataDsamp = zeros(14,numSampVA);
for i = 1:16;
    dataDsamp(i,:) = data(:,i);
end
dataVA = zeros(numApl,numSampVA/10);
s=[2 3 4 7 8 16];% 8 10 11 13 16];
for i = 1:6;
    %figure
    dataVA(i,1:numSampVA/10) = downsample(dataDsamp(s(i),:),10); 
    %plot(1:numSampVA/10,dataVA(i,1:numSampVA/10));
    %axis('tight')
    %pause
    %close all
end
%break
numSampVA=numSampVA/10;

%break
numState=24;                %total states of the 3 appliances
y = sum(dataVA(:,:),1);
numTrial=1;
dataEsVa = zeros(numApl,numSampVA,numTrial);
acc = zeros(2,numTrial);
ratingVA = [3; 9; 92; 410;
    70; 200;300; 400;500;600;700;
    4; 65;120;180;
    17;80;450;550;1620;
    1580;
    40;70;2700];
ac(1,:)=zeros(1,numState);  ac(1,4)=1; ac(1,10)=1; ac(1,14)=1;
ac(2,:)=zeros(1,numState);  ac(2,4)=1; ac(2,8)=1; ac(2,15)=1;
ac(3,:)=zeros(1,numState);  ac(3,3)=1; ac(3,19)=1;
ac(4,:)=zeros(1,numState);  ac(4,3)=1; ac(4,20)=1;% 
ac(5,:)=zeros(1,numState);  ac(5,9:11)=1;ac(5,21)=1;
ac(6,:)=zeros(1,numState); ac(6,4)=1; ac(6,14:15)=1;

h = ones(numApl,1);

F = zeros(numApl,numState);
F(1,1:4)=ones(1,4);
F(2,5:11)=ones(1,7);
F(3,12:15)=ones(1,4);
F(4,16:20)=ones(1,5);
F(5,21)=1;
F(6,22:24)=ones(1,3);
u1 = ones(5,1);

%% Integer LP initialization begins
c = [1; zeros(numState,1)];
b = [1; zeros(numState,1)];
Fh=F'*h;
a = [0; (ratingVA.*Fh)];
B1=[];
B2 = [0 ac(1,:); 0 ac(2,:); 0 ac(3,:); 0 ac(4,:); 
    0 ac(5,:);0 ac(6,:);];
Beq = [zeros(5,1) [F(1,:); F(2,:); F(3,:); F(4,:); F(6,:);]];
A1 = [-(a+b)'; (a-b)'; Beq];
A2 = [-(a+b)'; (a-b)'; B2];
intcon = 2:numState+1;
lb = zeros(numState+1,1); ub = [inf; ones(numState,1)];
cpu_IP = 0; cpu_ALIP=0;
for trial = 1:numTrial;
    f = ones(numState,numSampVA,2);
    zf = zeros(numApl,numSampVA,2);
    for i = 2:numSampVA;       
        options=optimoptions(@intlinprog,'display','off');
        tm = cputime;
        eb = [-y(1,i);  y(1,i); 1;1;1;1;1;];
        zhat = intlinprog(c,intcon,A1,eb,[],[],lb,ub,options);
        f(:,i,1)  = round(zhat(2:end,1));%
        zf(:,i,1) = F*(f(:,i,1).*ratingVA);
        elt = cputime-tm;
        cpu_IP = cpu_IP+elt;
        
        tm = cputime;
        eb = [-y(1,i);  y(1,i); 1;1;1;1;1;1;];
        zhat = intlinprog(c,intcon,A2,eb,Beq,u1,lb,ub,options);
        f(:,i,2)  = round(zhat(2:end,1));%
        zf(:,i,2) = F*(f(:,i,2).*ratingVA);
        if i>17 && zf(1,i,2)==410 
            if median(zf(1,i-17:i,2))==9;
                zf(1,i,2) = 9;
            elseif median(zf(1,i-17:i,2))==3;
                zf(1,i,2) = 3;
            end
        end
        if i>3 && zf(1,i,2)==92
            if median(zf(1,i-3:i,2))==9;
                zf(1,i,2) = 9;
            elseif median(zf(1,i-3:i,2))==3;
                zf(1,i,2) = 3;
            end
        end
        if zf(1,i-1,2)==410 && zf(1,i,2)<10;
            zf(1,i-1,2)=92;
        end
        if zf(2,i-1,2)==700 && zf(2,i,2)~=70;
            zf(2,i,2) = 700;
        end
        if zf(2,i-1,2)==500 && zf(2,i,2)~=70;
            zf(2,i,2) = 500;
        end
        if i>15 && zf(2,i-15,2)~=700 && median(zf(2,i-15:i,2))==700;
            zf(2,i-15,2) = 700;
        end
        if i>15 && zf(2,i-15,2)~=600 && median(zf(2,i-15:i,2))==600;
            zf(2,i-15,2) = 600;
        end
        if i>15 && zf(3,i-15,2)==180
            kk = median(zf(3,i-15:i,2));
            if kk==4;
                zf(3,i-15,2) = 4;
            end
            if kk==120;
                zf(3,i-15,2)=120;
            end
            if kk==65;
                zf(3,i-15,2)=65;
            end
        end
        if i>15 && zf(3,i-15,2)==120
            kk= median(zf(3,i-15:i,2));
            if kk==4;
                zf(3,i-15,2) = 4;
            end
            if kk==65;
                zf(3,i-15,2)=65;
            end
        end
        if zf(4,i-1,2)==550 && zf(4,i,2)==1620;
            zf(4,i,2)=550;
        end
        if zf(4,i-1,2)==450 && zf(4,i,2)==1620;
            zf(4,i,2)=450;
        end
        if i>5 && zf(5,i-5,2)==1580 && median(zf(5,i-5:i,2))==0;
            zf(5,i-5,2)=0;
        end
        if i>7 && zf(5,i-7,2)==0 && median(zf(5,i-7:i,2))==1580;
            zf(5,i-7,2)=1580;
        end
        elt = cputime-tm;
        cpu_ALIP = cpu_ALIP+elt;
    end
    for i=2:numSampVA
        tm = cputime;
        zft=[];
        zfn=[]; zfx=[];
        if zf(5,i,2)==1580;
            zft = 1580;
            zfn = 1580; 
            zfx = 1680;
        end
        if zf(6,i,2)==2700;
            zft = [zft 2700];
            zfn = [zfn 2700];
            zfx = [zfx 2800];            
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
            for k = 5:6;
                if zft(j)==zf(k,i,2);
                    zf(k,i,2) = zhat(j+1,1);
                end
            end
        end
        elt = cputime-tm;
        cpu_ALIP =cpu_ALIP+elt;
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
    for i=1:6;    
        Sacc(i,trial,1)=1-sum(abs(e1(i,:,trial)),1)/(2*sum(abs(dataVA(i,:,trial)),1));
        Sacc(i,trial,2)=1-sum(abs(e2(i,:,trial)),1)/(2*sum(abs(dataVA(i,:,trial)),1));  
    end
 end
 IP_AC = Sacc(:,1,1)
 ALIP_AC = Sacc(:,1,2)