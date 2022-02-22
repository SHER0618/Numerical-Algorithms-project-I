%search eigenvalue in a real schur
%input: a real schur matrix H
%output:A n*2 matrix D
function D=eig_search(H)
n=size(H,1);
D=zeros(n,2);
i=1;
while(i<n+1)
    if(i==n||H(i+1,i)==0)
        D(i,1)=H(i,i);
        i=i+1;
    else
        t=H(i,i)+H(i+1,i+1);
        s=H(i,i)*H(i+1,i+1)-H(i,i+1)*H(i+1,i);
        delta=t^2-4*s;
        if(delta>=0)
            D(i,1)=(t+sqrt(abs(delta)))/2;
            D(i+1,1)=(t-sqrt(abs(delta)))/2;
        else
            D(i,1)=t/2;
            D(i+1,1)=D(i,1);
            D(i+1,2)=sqrt(abs(delta))/2;
            D(i,2)=-D(i+1,2);
        end
        i=i+2;
    end
end
end