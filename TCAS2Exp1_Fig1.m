clear all
close all
load Exp1ResultSM.mat
        e1 =dataVA(:,:,1:72)-dataEsVa(:,:,1:72,1);
        e2 =dataVA(:,:,1:72)-dataEsVa(:,:,1:72,2);
        numTrial=72;
 for trial=1:numTrial;   
    for i=1:4;    
        Sacc(i,trial,1)=1-sum(abs(e1(i,:,trial)),1)/(2*sum(abs(dataVA(i,:,trial)),1));
        Sacc(i,trial,2)=1-sum(abs(e2(i,:,trial)),1)/(2*sum(abs(dataVA(i,:,trial)),1));
    end
 end
figure;
subplot(2,2,1);
plot(1:numTrial,Sacc(1,:,1),'-r','linewidth',2,'MarkerSize',10);
hold on
plot(1:numTrial,Sacc(1,:,2),'-b','linewidth',2,'MarkerSize',10);
%title('CDE')
%legend('IP', 'ALIP')
xlabel('Number of blocks');
ylabel('AC, CDE');
axis('tight')
grid on
subplot(2,2,2);
plot(1:numTrial,Sacc(2,:,1),'-r','linewidth',2,'MarkerSize',10);
hold on
plot(1:numTrial,Sacc(2,:,2),'-b','linewidth',2,'MarkerSize',10);
%title( 'FRG')
%legend('IP', 'ALIP')
xlabel('Number of blocks');
ylabel('AC, FRG');
axis('tight')
grid on
subplot(2,2,3);
plot(1:numTrial,Sacc(4,:,1),'-r','linewidth',2,'MarkerSize',10);
hold on
plot(1:numTrial,Sacc(4,:,2),'-b','linewidth',2,'MarkerSize',10)
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
plot(1:numTrial,Sacc(3,:,1),'-r','linewidth',1.5,'MarkerSize',10);
hold on
plot(1:numTrial,Sacc(3,:,2),'-b','linewidth',1.5,'MarkerSize',10);
title('HPE')
legend('IP', 'PIP')
xlabel('Number of blocks');
ylabel('AC');
axis('tight')
grid on
