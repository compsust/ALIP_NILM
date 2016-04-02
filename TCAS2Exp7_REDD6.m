close all
clear all
numSampVA = 192190;     % number of samples
Ni = 1;
rngS = [Ni,1,Ni-1+numSampVA,13]; % range
%rngS = [1,1,numSampVA,12]; % range
datafile = ['REDD\REDDhouse6_3sec_VA.csv'; ];

%% generating energy meter output y
numApl = 7; % number of appliances
data = csvread(datafile,Ni,1,rngS);
dataDsamp = zeros(13,numSampVA);
for i = 1:13;
    dataDsamp(i,:) = data(:,i);
end
dataVA = zeros(numApl,numSampVA/10);
s=[5 6 7 9 11 12 13];
for i = 1:7;
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
numState=20;                
y = sum(dataVA(:,:),1);
numTrial=1;
dataEsVa = zeros(numApl,numSampVA,numTrial);
acc = zeros(2,numTrial);
ratingVA = [530;710;
    3;950;
    145; 410;
    50;145;200;250;380;
    65;110;130;1000;
    400;2300;
    50;240;800;];
ac(1,:)=zeros(1,numState);  ac(1,4)=1; ac(1,7)=1; ac(1,12)=1;
ac(1,18)=1;
ac(2,:)=zeros(1,numState);  ac(2,1)=1; ac(2,6)=1; ac(2,7)=1; ac(2,12)=1;
ac(3,:)=zeros(1,numState);  ac(3,2)=1; ac(3,4)=1;ac(3,6)=1;
ac(3,11)=1;ac(3,15)=1;ac(3,17)=1;ac(3,20)=1;
ac(4,:)=zeros(1,numState);  ac(4,19)=1; ac(4,16)=1; 
ac(5,:)=zeros(1,numState); ac(5,19)=1;ac(5,5)=1; ac(5,8)=1;

h = ones(numApl,1);

F = zeros(numApl,numState);
F(1,1:2)=ones(1,2);
F(2,3:4)=ones(1,2);
F(3,5:6)=ones(1,2);
F(4,7:11)=ones(1,5);
F(5,12:15)=ones(1,4);
F(6,16:17)=ones(1,2);
F(7,18:20)=ones(1,3);
u1 = ones(3,1);

%% Integer LP initialization begins
c = [1; zeros(numState,1)];
b = [1; zeros(numState,1)];
Fh=F'*h;
a = [0; (ratingVA.*Fh)];
B1 = [0 F(1,:); 0 F(3,:); 0 F(4,:); 0 F(6,:);];
B2 = [0 F(1,:); 0 F(3,:); 0 F(4,:); 0 F(6,:); 0 ac(1,:);
    0 ac(2,:);0 ac(3,:); 0 ac(4,:);0 ac(5,:);];
Beq = [zeros(3,1) [F(2,:); F(5,:); F(7,:); ]];
A1 = [-(a+b)'; (a-b)';B1;Beq];
A2 = [-(a+b)'; (a-b)';B2];
intcon = 2:numState+1;
lb = zeros(numState+1,1); ub = [inf; ones(numState,1)];
cpu_IP=0; cpu_ALIP=0;
for trial = 1:numTrial;
    f = ones(numState,numSampVA,2);
    zf = zeros(numApl,numSampVA,2);
    for i = 2:numSampVA;
        options=optimoptions(@intlinprog,'display','off');
        tm = cputime;       
        eb = [-y(1,i);  y(1,i); 1;1;1;1;1;1;1;];
        zhat = intlinprog(c,intcon,A1,eb,[],[],lb,ub,options);
        f(:,i,1)  = round(zhat(2:end,1));%
        zf(:,i,1) = F*(f(:,i,1).*ratingVA);
        elt = cputime-tm;
        cpu_IP = cpu_IP+elt;
        
        tm = cputime;
        eb = [-y(1,i);  y(1,i); 1;1;1;1;1;1;1;1;1;];
        zhat = intlinprog(c,intcon,A2,eb,Beq,u1,lb,ub,options);
        f(:,i,2)  = round(zhat(2:end,1));%
        zf(:,i,2) = F*(f(:,i,2).*ratingVA);
        
        if i>9 && zf(1,i-9,2)==530 && median(zf(1,i-9:i,2))==0;
           zf(1,i-9,2) = 0;
        end
        if zf(3,i-1,2)==0 && zf(3,i,2)==410;
            zf(3,i,2)=145;
        end
        if zf(3,i-1,2)==410 && zf(3,i,2)==0;
            zf(3,i,2)=145;
        end
        if i>9 && zf(3,i-9,2)==0 && median(zf(3,i-9:i,2))==145;
           zf(3,i-9,2) = 145;
        end
        if i>9 && zf(4,i-9,2)==145 && median(zf(4,i-9:i),2)==0;
            zf(4,i-9,2) = 0;
        end
        if i>7 && zf(4,i-7,2)==200 && median(zf(4,i-7:i),2)==0;
            zf(4,i-7,2) = 0;
        end
        if i>7 && zf(4,i-7,2)==250 && median(zf(4,i-7:i),2)==0;
            zf(4,i-7,2) = 0;
        end
        if i>9 && zf(4,i-9,2)==380 && median(zf(4,i-9:i),2)==200;
            zf(4,i-9,2) = 200;
        end
        if i>9 &&  median(zf(4,i-9:i),2)==50;
            if zf(4,i-9,2)==250 
                zf(4,i-9,2) = 50;
            end
            if zf(4,i-9,2)==200
                zf(4,i-9,2) = 50;
            end
            if zf(4,i-9,2)==145;
                zf(4,i-9,2) = 50;
            end
        end
        if i>19 && zf(4,i-19,2)==0 && median(zf(4,i-19:i),2)==200;
            zf(4,i-19,2) = 200;
        end
        if i>19 && zf(7,i-19,2)==240 && median(zf(7,i-19:i),2)==50;
            zf(7,i-19,2) = 50;
        end
        elt = cputime-tm;
        cpu_ALIP = cpu_ALIP+elt;
    end
    for i=2:numSampVA
        tm = cputime;
        zft=[];
        zfn=[]; zfx=[];
        if zf(3,i,2)==145;
            zft = 145;
            zfn = 145; 
            zfx = 160;
        end
        if zf(6,i,2)==2300;
            zft = [zft 2300];
            zfn = [zfn 2300];
            zfx = [zfx 2350];            
        end
        if zf(6,i,2)==400;
            zft = [zft 400];
            zfn = [zfn 400];
            zfx = [zfx 450];            
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