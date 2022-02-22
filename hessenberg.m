%hessenberg a matrix A
%input:a matrix A
%output:a hessenberg matrix H=QTAQ,the ith household w exists in W ith
%column
function [Q,H]=hessenberg(H)
    [n,~]=size(H);
    Q=eye(n,n);
    for i=1:n-2
        w=household(H(i+1:n,i));
        H(i+1:n,i:n)=H(i+1:n,i:n)-2*w*(w'*H(i+1:n,i:n));
        H(1:n,i+1:n)=H(1:n,i+1:n)-2*(H(1:n,i+1:n)*w)*w';
        Q(1:n,i+1:n)=Q(1:n,i+1:n)-2*(Q(1:n,i+1:n)*w)*w';
        H(i+2:n,i)=0;
    end
end