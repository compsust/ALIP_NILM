clear all
close all
numSampVA = 406740;     % number of samples
Ni = 1;
rngS = [Ni,1,Ni-1+numSampVA,12]; % range
%rngS = [1,1,numSampVA,12]; % range
datafile = ['REDD\REDDhouse1_3sec_VA.csv'; ];

%% generating energy meter output y
numApl = 7; % number of appliances
data = csvread(datafile,Ni,1,rngS);
dataDsamp = zeros(12,numSampVA);
for i = 1:12;
    dataDsamp(i,:) = data(:,i);
end

numSampVA = numSampVA/20;
dataVA = zeros(numApl,numSampVA);
s=[2 3 4 8 7 9 12];%[2 3 4 6 7 8 9 12];
for i = 1:numApl;
    %figure
    dataVA(i,1:numSampVA) = downsample(dataDsamp(s(i),:),20); 
    %plot(1:numSampVA,dataVA(i,:));
    %axis('tight')
    %pause
    %close all
end
%break
numState=13;                %total states of the 3 appliances
numTrial=1;
dataEsVa = zeros(numApl,numSampVA,numTrial);
acc = zeros(2,numTrial);
y = sum(dataVA(:,:),1);

ratingVA = [4120; 
    6; 186; 425;
    1096; 215;
    1550;
    2740;550;
    1600;
    36; 80; 3200;];
ac(1,:)=[0 0 0 1 1 0 0 0 0 0 0 0 0];
ac(2,:)=[1 0 0 0 1 0 1 0 0 1 0 0 0];
h = ones(numApl,1);

%% current ratings of the states of the selected appliance is in r
F = zeros(numApl,numState);
F(1,1)=ones(1,1);
F(2,2:4)=ones(1,3);
F(3,5:6)=ones(1,2);
F(4,7)=ones(1,1);
F(5,8:9)=ones(1,2);
F(6,10)=ones(1,1);
F(7,11:13)=ones(1,3);

u1 = ones(2,1);

%% Integer LP initialization begins
c = [1; zeros(numState,1)];
b = [1; zeros(numState,1)];
Fh=F'*h;
a = [0; (ratingVA.*Fh)];
B1 = [zeros(1,1) F(3,:);zeros(1,1) F(5,:);];
B2 = [zeros(1,1) F(3,:);zeros(1,1) F(5,:);
     [zeros(1,1) ac(1,:)];[zeros(1,1) ac(2,:)];]; %
Beq = [zeros(2,1) [F(2,:); F(7,:)]];
A2 = [-(a+b)'; (a-b)'; B2];
A1 = [-(a+b)'; (a-b)'; B1; Beq];

intcon = 2:numState+1;
lb = zeros(numState+1,1); ub = [inf; ones(numState,1)];
cr = [1;zeros(numApl,1)];
br = [1;zeros(numApl,1)];
cpu_IP=0; cpu_ALIP=0;
for trial=1:numTrial
    f = zeros(numState,numSampVA,2);
    zf = zeros(numApl,numSampVA,2);
    for i = 2:numSampVA;
        options=optimoptions(@intlinprog,'display','off');
        tm = cputime;
        eb = [-y(1,i);  y(1,i); 1;1;1;1;];
        zhat = intlinprog(c,intcon,A1,eb,[],[],lb,ub,options);
        f(:,i,1)  = round(zhat(2:end,1));%
        zf(:,i,1) = F*(f(:,i,1).*ratingVA);
        elt = cputime-tm;
        cpu_IP= cpu_IP+elt;
        tm = cputime;
        eb = [-y(1,i);  y(1,i); 1;1;1;1];
        zhat = intlinprog(c,intcon,A2,eb,Beq,u1,lb,ub,options);
        f(:,i,2)  = round(zhat(2:end,1));%
        zf(:,i,2) = F*(f(:,i,2).*ratingVA);
        if zf(2,i-1,2)==6 && zf(2,i,2)==425;
            zf(2,i,2) = 186;
        end   
        if zf(2,i,2)==6 && zf(2,i-1,2)==425;
            zf(2,i,2) = 186;
        end
        if i>5 && zf(2,i-5,2)==425 && median(zf(2,i-5:i,2))<425;
            zf(2,i-5,2)=186;
        end
        if  zf(3,i-1,2)==0 && zf(3,i,2)==215 ;
            zf(3,i,2)=0;
        end
        if i>5 && zf(3,i-5,2)==0 && median(zf(3,i-5:i,2))<1096;
            zf(3,i-4,2)=0;
        end
        if i>5 && zf(5,i-5,2)==0 && median(zf(5,i-5:i,2))<550;
            zf(5,i-4,2)=0;
        end
        elt = cputime-tm;
        cpu_ALIP = cpu_ALIP+elt;
    end
    for i = 2:numSampVA;
        tm = cputime;
        zft=[];
        zfn=[]; zfx=[];
        if zf(2,i,2)==186;
            zft = 186;
            zfn = 186; 
            zfx = 200;
        end
        if zf(5,i,2)==2740;
            zft = [zft 2740];
            zfn = [zfn 2740];
            zfx = [zfx 3100];            
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
            for k = 1:5;
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