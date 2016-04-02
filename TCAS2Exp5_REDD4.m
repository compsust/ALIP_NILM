clear all
close all
numSampVA = 428070;     % number of samples
Ni = 1;
rngS = [Ni,1,Ni-1+numSampVA,13]; % range
%rngS = [1,1,numSampVA,12]; % range
datafile = ['REDD\REDDhouse4_3sec_VA.csv'; ];

%% generating energy meter output y
numApl = 8; % number of appliances
data = csvread(datafile,Ni,1,rngS);
dataDsamp = zeros(13,numSampVA);
for i = 1:13;
    dataDsamp(i,:) = data(:,i);
end
dataVA = zeros(numApl,numSampVA/10);
s=[2 4 6 7 8 10 11 12];
for i = 1:numApl;
    %figure
    dataVA(i,1:numSampVA/10) = downsample(dataDsamp(s(i),:),10); 
    %plot(1:numSampVA/10,dataVA(i,1:numSampVA/10));
    %axis('tight')
    %pause
    %close all
end
%break
numSampVA=numSampVA/10;

numState=20;                %total states of the 3 appliances
y = sum(dataVA(:,:),1);
numTrial=1;
dataEsVa = zeros(numApl,numSampVA,numTrial);
acc = zeros(2,numTrial);
ratingVA = [89; 174; 240; 308;
    5;10; 130;580;
    250;780;1016;
    5;400;1900;
    10;50;
    1;6;
    1320;
    1150];
ac(1,:)=zeros(1,numState); ac(1,4)=1;ac(1,9)=1; ac(1,16)=1; ac(1,18)=1;
ac(2,:)=zeros(1,numState); ac(2,2)=1; ac(2,7)=1; ac(2,16)=0;
ac(3,:)=zeros(1,numState); ac(3,1:4)=1;ac(3,20)=1;
ac(4,:)=zeros(1,numState); ac(4,10)=1; ac(4,4)=1;ac(4,13)=1;ac(4,16)=1;
ac(4,18)=1; ac(4,6)=1;
ac(5,:)=zeros(1,numState);ac(5,11)=1;ac(5,13)=1;ac(5,16)=1;
ac(6,:)=zeros(1,numState);ac(6,11)=1;ac(6,8)=1;ac(6,4)=1;
ac(7,:)=zeros(1,numState);ac(7,11)=1;ac(7,8)=1;ac(7,4)=1;
ac(8,:)=zeros(1,numState);ac(8,4)=1;ac(8,8)=1;ac(8,13)=1;
ac(8,16)=1;
ac(9,:)=zeros(1,numState);ac(9,20)=1;ac(9,11)=1;ac(9,8)=1;
ac(9,4)=1;ac(9,15)=1;%
ac(10,:)=zeros(1,numState);ac(10,7)=1;ac(10,16)=0;ac(10,2)=1;
ac(11,:)=zeros(1,numState);ac(11,1)=1;ac(11,16)=1;
ac(12,:)=zeros(1,numState); ac(12,5)=1; ac(12,18)=1;
ac(13,:)=zeros(1,numState); ac(13,11)=1; ac(13,14)=1;

h = ones(numApl,1);

F = zeros(numApl,numState);
F(1,1:4)=ones(1,4);
F(2,5:8)=ones(1,4);
F(3,9:11)=ones(1,3);
F(4,12:14)=ones(1,3);
F(5,15:16)=ones(1,2);
F(6,17:18)=ones(1,2);
F(7,19)=ones(1,1);
F(8,20)=ones(1,1);
u1 = ones(2,1);

%% Integer LP initialization begins
c = [1; zeros(numState,1)];
b = [1; zeros(numState,1)];
Fh=F'*h;
a = [0; (ratingVA.*Fh)];
B1 = [0 F(1,:); 0 F(3,:);0 F(5,:); 0 F(6,:);];
B2 = [0 F(1,:); 0 F(3,:);0 F(5,:); 0 F(6,:); 0 ac(1,:); 0 ac(2,:);
    0 ac(3,:); 0 ac(4,:);0 ac(5,:); 0 ac(6,:);0 ac(7,:);0 ac(8,:);
    0 ac(9,:);0 ac(10,:);0 ac(11,:);0 ac(12,:);0 ac(13,:);];
Beq = [zeros(2,1) [F(2,:); F(4,:);]];
A1 = [-(a+b)'; (a-b)'; B1; Beq];
A2 = [-(a+b)'; (a-b)'; B2];
intcon = 2:numState+1;
lb = zeros(numState+1,1); ub = [inf; ones(numState,1)];
cpu_IP=0; cpu_ALIP=0;
for trial = 1:numTrial;
    f = ones(numState,numSampVA,2);
    zf = zeros(numApl,numSampVA,2);
    for i = 2:numSampVA;
        options=optimoptions(@intlinprog,'display','off');
        tm=cputime;
        eb = [-y(1,i);  y(1,i); 1;1;1;1;1;1;];
        zhat = intlinprog(c,intcon,A1,eb,[],[],lb,ub,options);
        f(:,i,1)  = round(zhat(2:end,1));%
        zf(:,i,1) = F*(f(:,i,1).*ratingVA);
        elt = cputime-tm;
        cpu_IP = cpu_IP+elt;
        
        tm = cputime;
        eb = [-y(1,i);  y(1,i); 1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;];
        zhat = intlinprog(c,intcon,A2,eb,Beq,u1,lb,ub,options);
        f(:,i,2)  = round(zhat(2:end,1));%
        zf(:,i,2) = F*(f(:,i,2).*ratingVA);
        if i>5 && zf(1,i-5,2)==89 && median(zf(1,i-5:i,2))==0;
            zf(1,i-5,2) = 0;
        end
        if i>63 && zf(1,i-63,2)==174 && median(zf(1,i-63:i,2))==0;
            zf(1,i-63,2) = 0;
        end
        if i>37 && zf(1,i-37,2)==240 && median(zf(1,i-37:i,2))==0;
            zf(1,i-37,2) = 0;
        end
        if i>7 && zf(1,i-7,2)==308 && median(zf(1,i-7:i,2))==0;
            zf(1,i-7,2) = 0;
        end
        if i>25 && zf(1,i-25,2)==308 && median(zf(1,i-25:i,2))==240;
            zf(1,i-25,2) = 89;
        end
        if i>25 && zf(1,i-25,2)==308 && median(zf(1,i-25:i,2))==89;
            zf(1,i-25,2) = 89;
        end
        if i>25&& zf(1,i-25,2)==240 && median(zf(1,i-25:i,2))==89;
            zf(1,i-25,2) = 89;
        end
        if zf(2,i-1,2)==5 && zf(2,i,2)==580;
            zf(2,i,2) = 5;
        end
        if zf(2,i-1,2)==10 && zf(2,i,2)==580;
            zf(2,i,2) = 10;
        end
        if zf(2,i-1,2)==580 && zf(2,i,2)==130;
            zf(2,i-1,2) = 130;
        end
        if i>5 && zf(2,i-5,2)==5 && median(zf(2,i-5:i,2))==130;
            zf(2,i-5,2)=130;
        end
        if i>25 && zf(3,i-25,2)==250 && median(zf(3,i-25:i,2))==0;
            zf(3,i-25,2) = 0;
        end
        if i>5 && zf(3,i-5,2)==1016 && median(zf(3,i-5:i,2))==0;
            zf(3,i-5,2) = 0;
        end
        if i>11 && zf(4,i-11,2)==400 && median(zf(4,i-11:i,2))==5;
            zf(4,i-11,2) = 5;
        end
        if i>45 && zf(5,i-45,2)==10 && median(zf(5,i-45:i,2))==0;
            zf(5,i-45,2) = 0;
        end
        if i>45 && zf(5,i-45,2)==10 && median(zf(5,i-45:i,2))==10;
            zf(5,i-45,2) = 0;
        end
        if i>35 && zf(5,i-35,2)==50 && median(zf(5,i-35:i,2))==50;
            zf(5,i-35,2) = 0;
        end
        if i>45;
            mzf6 = median(zf(6,i-45:i,2));
        end
        if i>45 && zf(6,i-45,2)==6 && mzf6==6;
            zf(6,i-45,2) = 0;
        end
        if i>45 && zf(6,i-45,2)==1 && mzf6==0;
            zf(6,i-45,2) = 0;
        end
        if i>45 && mzf6==1 && zf(6,i,2)==6;
            zf(6,i,2)=1;
        end
        if i>45 && mzf6==0 && zf(6,i,2)==6;
            zf(6,i,2)=0;
        end
        elt = cputime-tm;
        cpu_ALIP = cpu_ALIP+elt;
    end
    for i=2:numSampVA
        tm = cputime;
        zft=[];
        zfn=[]; zfx=[];
        if zf(1,i,2)==89;
            zft = 89;
            zfn = 89; 
            zfx = 103;
        end
        if zf(1,i,2)==174;
            zft = [zft 174];
            zfn = [zfn 174];
            zfx = [zfx 185];            
        end
        if zf(1,i,2)==240;
            zft = [zft 240];
            zfn = [zfn 240];
            zfx = [zfx 255];            
        end
        if zf(1,i,2)==308;
            zft = [zft 308];
            zfn = [zfn 275];
            zfx = [zfx 308];            
        end
        if zf(2,i,2)==130;
            zft = [zft 130];
            zfn = [zfn 130];
            zfx = [zfx 160];            
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
            for k = 1:2;
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
    for i=1:8;    
        Sacc(i,trial,1)=1-sum(abs(e1(i,:,trial)),1)/(2*sum(abs(dataVA(i,:,trial)),1));
        Sacc(i,trial,2)=1-sum(abs(e2(i,:,trial)),1)/(2*sum(abs(dataVA(i,:,trial)),1));  
    end
 end
 IP_AC = Sacc(:,1,1)
 ALIP_AC = Sacc(:,1,2)