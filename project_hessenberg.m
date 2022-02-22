%project_hessenberg 
time_my=zeros(20,1);
time_matlab=zeros(20,1);
for n=100:100:2000
    temp_time_my=zeros(5,1);
    temp_time_matlab=zeros(5,1);
    for lp=1:5
        A=rand(n,n);
        tic;
        [Q,H]=hessenberg(A,1);
        temp_time_my(lp,1)=toc;
        tic;
        [Q,H]=hess(A);
        temp_time_matlab(lp,1)=toc;
    end
    time_my((n-100)/100+1,1)=sum(temp_time_my)/5;
    time_matlab((n-100)/100+1,1)=sum(temp_time_matlab)/5;
end
plot([1:20],time_my);
hold on;
plot([1:20],time_matlab);
title('The running time of hessenberg');
xlabel('The order of matrix');
ylabel('The runing time');
legend('My hessenberg','Matlab hessenberg');