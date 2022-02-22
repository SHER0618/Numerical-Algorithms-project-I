%choose eigenvalue
%input:E is a sorted eigenvalue from smalll to largest, k is the number we
%needed and k<n
%output:a E1 with real eigenvalue and a E2 with complex eigenvalue
function [E1,E2,t,m]=choose(E)
n=size(E,1);
E1=zeros(n,1);
E2=zeros(n,1);
t=1;%t -E1
m=1;%m -E2
for i=1:n
    if(imag(E(n-i+1,1))==0)
        E1(t,1)=real(E(n-i+1,1));
        t=t+1;
    else
        E2(m,1)=E(n-i+1,1);
        m=m+1;
    end
end
E1=E1(1:t-1,1);
E2=E2(1:m-1,1);
t=t-1;
m=m-1;