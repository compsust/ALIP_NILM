clear all
close all
numSampVA = 316840;     % number of samples
Ni = 1;
rngS = [Ni,1,Ni-1+numSampVA,10]; % range
%rngS = [1,1,numSampVA,12]; % range
datafile = ['REDD\REDDhouse2_3sec_VA.csv'; ];

%% generating energy meter output y
numApl = 6; % number of appliances
data = csvread(datafile,Ni,1,rngS);
dataDsamp = zeros(10,numSampVA);
for i = 1:10;
    dataDsamp(i,:) = data(:,i);
end
dataVA = zeros(numApl,numSampVA/5);
s=[2 3 8 5 4 6];
for i = 1:numApl;
    %figure
    dataVA(i,1:numSampVA/5) = downsample(dataDsamp(s(i),:),5); 
    %plot(1:numSampVA/10,dataVA(i,1:numSampVA/10));
    %axis('tight')
    %pause
    %close all
end
%break
numSampVA=numSampVA/5;

numState=17;               %total states of the 4 appliances
numTrial=1;
dataEsVa = zeros(numApl,numSampVA,numTrial);
acc = zeros(2,numTrial);
y = sum(dataVA(:,:),1);

ratingVA = [1050; 770;14; 
    8; 20; 42; 75; 110;150;
    2;250; 1200;
    5; 47; 1800;
    410;
    2];
ac(1,:) =[1 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0];
ac(2,:) =[1 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0];
ac(3,:) =[0 0 0 0 1 0 0 0 0 0 0 1 0 0 0 0 0];
ac(4,:) =[0 0 0 0 0 1 0 0 0 0 0 1 0 0 0 0 0];
ac(5,:) =[0 0 0 0 0 0 1 0 0 0 0 1 0 0 0 0 0];
ac(6,:) =[1 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0];
ac(7,:) =[1 0 0 0 0 0 0 0 1 0 0 0 0 1 0 0 0];
ac(8,:) =[0 0 0 0 0 1 0 0 0 0 0 0 0 1 0 0 0];
ac(9,:) =[0 0 0 0 1 0 0 0 0 0 0 0 0 1 0 0 0];
ac(10,:)=[0 0 0 0 0 0 0 1 0 0 0 0 0 1 0 0 0];
ac(11,:)=[0 1 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0];
ac(12,:)=[0 1 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0];
ac(13,:)=[0 0 0 0 0 0 0 0 1 0 0 0 0 0 1 0 0];
ac(14,:)=[1 0 0 0 0 0 0 1 1 0 1 0 0 0 0 1 0];%
ac(15,:)=[0 0 0 0 0 1 0 0 0 0 0 0 1 0 0 0 0];
ac(16,:)=[0 0 1 0 0 0 0 0 0 1 0 0 0 0 0 0 0];
ac(17,:)=[0 0 0 0 0 0 1 0 0 0 0 0 0 1 0 0 0];
ac(18,:)=[0 0 1 0 0 0 1 0 0 0 0 0 1 0 0 0 0];

h = ones(numApl,1);

%% current ratings of the states of the selected appliance is in r
F = zeros(numApl,numState);
F(1,1:3)=ones(1,3);
F(2,4:9)=ones(1,6);
F(3,10:12)=ones(1,3);
F(4,13:15)=ones(1,3);
F(5,16)=ones(1,1);
F(6,17)=ones(1,1);
u1 = ones(3,1);

%% Integer LP initialization begins
c = [1; zeros(numState,1)];
b = [1; zeros(numState,1)];
Fh=F'*h;
a = [0; (ratingVA.*Fh)];
B1 = [zeros(1,1) F(1,:); zeros(1,1) F(3,:);];
B2 = [zeros(1,1) F(1,:); zeros(1,1) F(3,:); 0 ac(1,:);
    0 ac(2,:); 0 ac(3,:); 0 ac(4,:); 0 ac(5,:);
    0 ac(6,:); 0 ac(7,:); 0 ac(8,:); 0 ac(9,:); 0 ac(10,:);
    0 ac(11,:); 0 ac(12,:); 0 ac(13,:); 0 ac(14,:);
    0 ac(15,:); 0 ac(16,:);0 ac(17,:); 0 ac(18,:)];
Beq = [zeros(3,1) [F(2,:); F(4,:); F(6,:);]];
A2 = [-(a+b)'; (a-b)'; B2];
A1 = [-(a+b)'; (a-b)'; B1; Beq];
intcon = 2:numState+1;
lb = zeros(numState+1,1); ub = [inf; ones(numState,1)];
cpu_IP=0;
cpu_ALIP=0;
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
        eb = [-y(1,i);  y(1,i); 1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;2];
        zhat = intlinprog(c,intcon,A2,eb,Beq,u1,lb,ub,options);
        f(:,i,2)  = round(zhat(2:end,1));%
        zf(:,i,2) = F*(f(:,i,2).*ratingVA);
        if i>5 && zf(1,i-5,2)==0 && median(zf(1,i-5:i,2))==14; 
            zf(1,i-5,2)=14;
        end
        if i>2 && zf(1,i-2,2)==770
            if median(zf(1,i-2:i,2))<14 ;
                zf(1,i-2,2)=0;
            elseif median(zf(1,i-2:i,2))<770;
                zf(1,i-2,2)=14;
            end
        end
        if i>5
            if zf(2,i-5,2)==20 && median(zf(2,i-5:i,2))==8;
                zf(2,i-5,2)= 8;
            end
        end
        if i>3
            if zf(2,i-3,2)==110 && median(zf(2,i-3:i,2))==8;
                zf(2,i-3,2)= 8;
            end
        end
        if i>9
            if zf(2,i-9,2)==8 && median(zf(2,i-9:i),2)==150;
                zf(2,i-9,2)= 150;
            end
        end
        if i>3 && zf(3,i-3,2)==1200
            if median(zf(3,i-3:i,2))<5; 
                zf(3,i-3,2)=0;
            elseif median(zf(3,i-3:i,2))<1200; 
                zf(3,i-3,2)=250;        
            end
        end
        if i>20 && (zf(3,i-20,2))==250 && median(zf(3,i-20:i,2))<5;
            zf(3,i-20,2)=0;
        end
        if i>5 && zf(5,i-5,2)==410 && median(zf(5,i-5:i,2)) == 0 ;
            zf(5,i-5,2) = 0;
        end
        elt = cputime-tm;
        cpu_ALIP=cpu_ALIP+elt;
    end
    for i=2:numSampVA
        tm = cputime;
        zft=[];
        zfn=[]; zfx=[];
        if zf(1,i,2)==1050;
            zft = 1050;
            zfn = 1050; 
            zfx = 1070;
        end
        if zf(2,i,2)==110;
            zft = [zft 110];
            zfn = [zfn 110];
            zfx = [zfx 115];            
        end
        if zf(2,i,2)==150;
            zft = [zft 150];
            zfn = [zfn 150];
            zfx = [zfx 170];            
        end
        if zf(3,i,2)==1200;
            zft = [zft 1200];
            zfn = [zfn 1185];
            zfx = [zfx 1230];            
        end
        if zf(4,i,2)==1800;
            zft = [zft 1800];
            zfn = [zfn 1800];
            zfx = [zfx 1850];            
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
            for k = 1:4;
                if zft(j)==zf(k,i,2);
                    zf(k,i,2) = zhat(j+1,1);
                end
            end
        end
        elt = cputime-tm;
        cpu_ALIP=cpu_ALIP+elt;
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