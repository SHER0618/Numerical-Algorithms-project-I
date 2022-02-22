%project-AED-running-time1
clear;
clc;
time_double=zeros(30,1);
time_sextuple=zeros(30,1);
time_AED=zeros(30,1);
for i=620:20:1200
        A=rand(i,i);
        tic;
        [E,H,Q]=double_shift_QR_algorithm(A,1);
        time_double((i-600)/20,1)=toc;
        tic;
        [E,H,Q]=sextuple_shift_QR_algorithm(A,1);
        time_sextuple((i-600)/20,1)=toc;
        tic;
        [E,Q,H]=aggressive_early_deflation(A,1);
        time_AED((i-600)/20,1)=toc;
end
index=[620:20:1200]';
plot(index,time_double);
hold on;
plot(index,time_sextuple);
hold on;
plot(index,time_AED);
title('The running time of traditional QR Algorithm and Aggressive Early Deflation');
xlabel('The order of matrix');
ylabel('The runing time');
legend('Double-Shift-QR-Algorithm','Sextuple-Shift-QR-Algorithm','Aggressive Early Deflation');