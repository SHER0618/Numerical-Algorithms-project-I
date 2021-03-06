%project_ill-conditioned matrix
A=-1+2*rand(500,500,23);
matrix_norm=zeros(23,1);
condition_number=zeros(23,1);
for i=1:23
    [S,V,D]=svd(A(1:500,1:500,i));
    V(1,1)=10^((i+5)/4);
    V(500,500)=10^((-i-5)/4);
    A(1:500,1:500,i)=S'*V*D;
    matrix_norm(i,1)=V(1,1);
    condition_number(i,1)=V(1,1)^2;
end
double_=zeros(23,1);
sextuple_=zeros(23,1);
AED_=zeros(23,1);
Matlab_=zeros(23,1);
double_Q=zeros(23,1);
sextuple_Q=zeros(23,1);
AED_Q=zeros(23,1);
Matlab_Q=zeros(23,1);
for i=1:23
    T=A(1:500,1:500,i);
    [E,H,Q]=double_shift_QR_algorithm(T,1);
    double_(i,1)=norm(T*Q-Q*H,'fro');
    double_Q(i,1)=norm(Q'*Q-eye(500),'fro');
    [E,H,Q]=sextuple_shift_QR_algorithm(T,1);
    sextuple_(i,1)=norm(T*Q-Q*H,'fro');
    sextuple_Q(i,1)=norm(Q'*Q-eye(500),'fro');
    [E,Q,H]=aggressive_early_deflation(T,1);
    AED_(i,1)=norm(T*Q-Q*H,'fro');
    AED_Q(i,1)=norm(Q'*Q-eye(500),'fro');
    [Q,H]=schur(T);
    Matlab_(i,1)=norm(T*Q-Q*H,'fro');
    Matlab_Q(i,1)=norm(Q'*Q-eye(500),'fro');
end
double_=double_./matrix_norm;
sextuple_=sextuple_./matrix_norm;
AED_=AED_./matrix_norm;
Matlab_=Matlab_./matrix_norm;
double_Q=double_Q./sqrt(500);
sextuple_Q=sextuple_Q./sqrt(500);
AED_Q=AED_Q./sqrt(500);
Matlab_Q=Matlab_Q./sqrt(500);
hold off;
loglog(condition_number,double_);
hold on;
loglog(condition_number,sextuple_);
hold on;
loglog(condition_number,AED_);
hold on;
loglog(condition_number,Matlab_);
xlabel('${\kappa _2}(A)$','Interpreter','latex','FontSize',18);
ylabel('$\frac{{{{\left\| {AQ - QH} \right\|}_F}}}{{{{\left\| A \right\|}_F}}}$','Interpreter','latex','rotation',0,'FontSize',18);
title('The numerical stability for schur decomposition','FontSize',12);
legend('Double-Shift-QR-Algorithm','Sextuple-Shift-QR-Algorithm','Aggressive Early Deflation','Matlab-schur');
hold off;
loglog(condition_number,double_Q);
hold on;
loglog(condition_number,sextuple_Q);
hold on;
loglog(condition_number,AED_Q);
hold on;
loglog(condition_number,Matlab_Q);
hold on;
legend('Double-Shift-QR-Algorithm','Sextuple-Shift-QR-Algorithm','Aggressive Early Deflation','Matlab-schur');
xlabel('${\kappa _2}(A)$','Interpreter','latex','FontSize',18);
ylabel('$\frac{{{{\left\| {{Q^T}Q - I} \right\|}_F}}}{{\sqrt n }}$','Interpreter','latex','rotation',0,'FontSize',15);
title('The numerical stability for orthogonal matrix','FontSize',12);