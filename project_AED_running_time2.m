%project-AED-running-time2
clear;
clc;
time_AED=zeros(30,1);
time_Matlab=zeros(30,1);
for i=1000:150:5350
    A=rand(i,i);
    tic;
    [E,Q,H]=aggressive_early_deflation(A,1);
    time_AED((i-850)/150,1)=toc;
    tic;
    [Q,H]=schur(A);
    time_Matlab((i-850)/150,1)=toc;
end
index=[1000:150:5350]';
plot(index,time_AED);
hold on;
plot(index,time_Matlab);
title('The running time of the Aggressive Early Deflation');
xlabel('The order of matrix');
ylabel('The runing time');
legend('Aggressive Early Deflation','Matlab-schur');