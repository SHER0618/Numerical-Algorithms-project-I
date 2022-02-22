%chasing QR
function [W0,H,m,k,Ava,least]=double_shift_chasing_iteration(H)
[n,~]=size(H);
m=floor(sqrt(n));
k=3*m;
%some_eigenvalue=size(D,1);
%and_eigenvalue=3*m-some_eigenvalue;
Ava=floor((n-3*m-1)/k);
W0=zeros(3*m+k+1,3*m+k+1,Ava+3);
%start
W=eye(3*m+1);
%E0=D(1:some_eigenvalue,1)+1i*D(1:some_eigenvalue,2);
%all_eigenvalue=eig(H(n-3*m+1:n,n-3*m+1:n));
[Eig,~,~]=aggressive_early_deflation(H(n-3*m+1:n,n-3*m+1:n),0);
%[E1,E2,~,~]=choose(eig(H(n-3*m+1:n,n-3*m+1:n)));
[E1,E2,~,~]=choose(Eig);
Eigenvalue=[E2;E1];
for i=1:m
     s=real(Eigenvalue(2*i-1,1)+Eigenvalue(2*i,1));
     t=real(Eigenvalue(2*i-1,1)*Eigenvalue(2*i,1));
     x=H(1,1)*H(1,1)+H(1,2)*H(2,1)-s*H(1,1)+t;
     y=H(2,1)*(H(1,1)+H(2,2)-s);
     z=H(2,1)*H(3,2);
     vector=[x;y;z];
     for Diana=0:3*(m-i)
         w=household(vector);
         q=max([1,Diana]);
         r=min([Diana+4,3*m+1]);
         H(Diana+1:Diana+3,q:(3*m+1))=H(Diana+1:Diana+3,q:(3*m+1))-2*w*(w'*H(Diana+1:Diana+3,q:(3*m+1)));
         H(1:r,Diana+1:Diana+3)=H(1:r,Diana+1:Diana+3)-2*(H(1:r,Diana+1:Diana+3)*w)*w';
         vector=H(Diana+2:Diana+4,Diana+1);
         W(1:3*m+1,Diana+1:Diana+3)=W(1:3*m+1,Diana+1:Diana+3)-2*(W(1:3*m+1,Diana+1:Diana+3)*w)*w';
     end
end
H(1:3*m+1,(3*m+2):n)=W'*H(1:3*m+1,(3*m+2):n);
W0(1:3*m+1,1:3*m+1,1)=W;
%middle
for i=1:Ava
    W=eye(3*m+k+1);  
    for t=1:k
        VEC=zeros(3,m);
        for pen=m:-1:1
            position=1+k*(i-1)+t-1+3*(pen-1);
            vector=H(position+1:position+3,position);
            w=household(vector);
            VEC(1:3,pen)=w;
            H(position+1:position+3,position:3*m+k*i+1)= H(position+1:position+3,position:3*m+k*i+1)- 2*w*(w'*H(position+1:position+3,position:3*m+k*i+1));
            H(1+k*(i-1):position+4,position+1:position+3)=H(1+k*(i-1):position+4,position+1:position+3)-2*(H(1+k*(i-1):position+4,position+1:position+3)*w)*w';
        end
        for plp=1:m
            w=VEC(1:3,plp);
            W(1:3*m+k+1,t+3*plp-2:t+3*plp)=W(1:3*m+k+1,t+3*plp-2:t+3*plp)-2*(W(1:3*m+k+1,t+3*plp-2:t+3*plp)*w)*w';
        end
    end
    H(1:k*(i-1),1+k*(i-1):3*m+k*i+1)=H(1:k*(i-1),1+k*(i-1):3*m+k*i+1)*W;
    H(1+k*(i-1):3*m+k*i+1,3*m+k*i+2:n)=W'*H(1+k*(i-1):3*m+k*i+1,3*m+k*i+2:n);
    W0(1:1+3*m+k,1:1+3*m+k,i+1)=W;
end
least=n-k*Ava-3*m-1;
W=eye(3*m+least+1);
for t=1:least
    VEC=zeros(3,m);
        for pen=m:-1:1
            position=1+k*i+t-1+3*(pen-1);
            vector=H(position+1:position+3,position);
            w=household(vector);
            VEC(1:3,pen)=w;
            H(position+1:position+3,position:n)= H(position+1:position+3,position:n)- 2*w*(w'*H(position+1:position+3,position:n));
            H(1+k*i:position+4,position+1:position+3)=H(1+k*i:position+4,position+1:position+3)-2*(H(1+k*i:position+4,position+1:position+3)*w)*w';
        end
        for plp=1:m
            w=VEC(1:3,plp);
            W(1:3*m+least+1,t+3*plp-2:t+3*plp)=W(1:3*m+least+1,t+3*plp-2:t+3*plp)-2*(W(1:3*m+least+1,t+3*plp-2:t+3*plp)*w)*w';
        end
end
H(1:k*i,1+k*i:n)=H(1:k*i,1+k*i:n)*W;
W0(1:1+3*m+least,1:1+3*m+least,Ava+2)=W;
%end:
W=eye(3*m+1);
for bella=m:-1:1
    position=n-3*m+3*(bella-1);
    for Carol=position:n-3
        r=min([Carol+4,n]);
        w=household(H(Carol+1:Carol+3,Carol));
        H(Carol+1:Carol+3,Carol:n)=H(Carol+1:Carol+3,Carol:n)-2*w*(w'*H(Carol+1:Carol+3,Carol:n));
        H(n-3*m:r,Carol+1:Carol+3)=H(n-3*m:r,Carol+1:Carol+3)-2*(H(n-3*m:r,Carol+1:Carol+3)*w)*w';
        W(1:3*m+1,Carol+2-n+3*m:Carol+4-n+3*m)=W(1:3*m+1,Carol+2-n+3*m:Carol+4-n+3*m)-2*(W(1:3*m+1,Carol+2-n+3*m:Carol+4-n+3*m)*w)*w';
    end
    w=household(H(n-1:n,n-2));
    H(n-3*m:n,n-1:n)=H(n-3*m:n,n-1:n)-2*(H(n-3*m:n,n-1:n)*w)*w';
    H(n-1:n,n-2:n)=H(n-1:n,n-2:n)-2*w*(w'*H(n-1:n,n-2:n));
    W(1:3*m+1,3*m:3*m+1)=W(1:3*m+1,3*m:3*m+1)-2*(W(1:3*m+1,3*m:3*m+1)*w)*w';
end
H(1:n-3*m-1,n-3*m:n)=H(1:n-3*m-1,n-3*m:n)*W;
W0(1:1+3*m,1:1+3*m,Ava+3)=W;
m=floor(sqrt(n));
k=3*m;
Ava=floor((n-3*m-1)/k);