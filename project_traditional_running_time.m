%project-traditional-running-time
clear;
clc;
time_double=zeros(90,1);
time_sextuple=zeros(90,1);
for i=110:10:1000
        A=rand(i,i);
        tic;
        [E,H,Q]=double_shift_QR_algorithm(A,1);
        time_double((i-100)/10,1)=toc;
        tic;
        [E,H,Q]=sextuple_shift_QR_algorithm(A,1);
        time_sextuple((i-100)/10,1)=toc;
end
index=[110:10:1000]';
plot(index,time_double);
hold on;
plot(index,time_sextuple);
title('The running time of traditional QR Algorithm');
xlabel('The order of matrix');
ylabel('The runing time');
legend('Double-Shift-QR-Algorithm','Sextuple-Shift-QR-Algorithm');