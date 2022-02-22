%double_shift_QR_algorithm
%input:A dense matrix H
%output:A n*2 matrix with all of the eigenvalue of A;
%The first column is the real part and the second column is the second
%column is the imaginary part, which the a+bi is front of a-bi and sort by
%the length of eigenvalues.
function [E,H,Q]=double_shift_QR_algorithm(A,flag)
    [n,~]=size(A);
    D=zeros(n,2);
    E=zeros(n,1);
    H=A;
    if(n==1)
        Q=eye(1);
        E(1,1)=A(1,1);
        return;
    end
    Q=eye(n);
    if flag==1
        %[Q,H]=hessenberg(H);
        [Q,H]=hess(H);
    end
    i=1;
    m=n;
    tol=1e-15;
    while(m-i+1>2)
        [W,H(i:m,i:m)]=double_shift_QR_iteration(H(i:m,i:m));
        [~,myu]=size(W);
        for tp=1:myu-1
            w=W(1:3,tp);
            Q(1:n,i-1+tp:i+1+tp)=Q(1:n,i-1+tp:i+1+tp)-2*(Q(1:n,i-1+tp:i+1+tp)*w)*w';
            H(1:i-1,i-1+tp:i+1+tp)=H(1:i-1,i-1+tp:i+1+tp)-2*(H(1:i-1,i-1+tp:i+1+tp)*w)*w';
            H(i-1+tp:i+1+tp,m+1:n)=H(i-1+tp:i+1+tp,m+1:n)-2*w*(w'*H(i-1+tp:i+1+tp,m+1:n));
        end
        w=W(1:2,myu);
        Q(1:n,i-1+myu:i+myu)=Q(1:n,i-1+myu:i+myu)-2*(Q(1:n,i-1+myu:i+myu)*w)*w';
        H(1:i-1,i-1+myu:i+myu)=H(1:i-1,i-1+myu:i+myu)-2*(H(1:i-1,i-1+myu:i+myu)*w)*w';
        H(i-1+myu:i+myu,m+1:n)=H(i-1+myu:i+myu,m+1:n)-2*w*(w'*H(i-1+myu:i+myu,m+1:n));
        while(m-i+1>2)
            py=0;
            if(abs(H(m,m-1))<tol)
                py=1;
            end
            if(abs(H(m-1,m-2))<tol)
                py=1;
            end
            if(abs(H(i+1,i))<tol)
                py=1;
            end
            if(abs(H(i+2,i+1))<tol)
                py=1;
            end
            if(py==0)
                break;
            end
            if abs(H(m,m-1))<tol
                H(m,m-1)=0;
                D(m,1)=H(m,m);
                m=m-1;
            else
                if abs(H(m-1,m-2))<tol
                    H(m-1,m-2)=0;
                    t=H(m,m)+H(m-1,m-1);
                    s=H(m,m)*H(m-1,m-1)-H(m,m-1)*H(m-1,m);
                    delta=t^2-4*s;
                    if(delta>=0)
                        D(m,1)=(t+sqrt(abs(delta)))/2;
                        D(m-1,1)=(t-sqrt(abs(delta)))/2;
                    else
                        D(m,1)=t/2;
                        D(m-1,1)=D(m,1);
                        D(m,2)=sqrt(abs(delta))/2;
                        D(m-1,2)=-D(m,2);
                    end
                    [temp_Q,H(m-1:m,m-1:m)]=schur(H(m-1:m,m-1:m));
                    Q(1:n,m-1:m)=Q(1:n,m-1:m)*temp_Q;
                    H(1:m-2,m-1:m)=H(1:m-2,m-1:m)*temp_Q;
                    H(m-1:m,m+1:n)=temp_Q'*H(m-1:m,m+1:n);
                    m=m-2;
                end
            end
            if abs(H(i+1,i))<tol
                H(i+1,i)=0;
                D(i,1)=H(i,i);
                i=i+1;
            else
                if abs(H(i+2,i+1))<tol
                    H(i+2,i+1)=0;
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
                    [temp_Q,H(i:i+1,i:i+1)]=schur(H(i:i+1,i:i+1));
                    Q(1:n,i:i+1)=Q(1:n,i:i+1)*temp_Q;
                    H(1:i-1,i:i+1)=H(1:i-1,i:i+1)*temp_Q;
                    H(i:i+1,i+2:n)=temp_Q'*H(i:i+1,i+2:n);
                    i=i+2;
                end
            end
        end
    end
%    idea=0;
 %   if(D(i,1)==0)
  %     idea=1;
%    end
 %   if(D(i+1,1)==0)
  %      idea=1;
%    end
 %   if(D(m,1)==0)
  %     idea=1;
   % end
    if(m>i)
        P=H(i:m,i:m);
        temp=eig(P);
        [U,H(i:m,i:m)]=schur(H(i:m,i:m));
        Q(1:n,i:m)=Q(1:n,i:m)*U;
        H(1:i-1,i:m)=H(1:i-1,i:m)*U;
        H(i:m,m+1:n)=U'*H(i:m,m+1:n);
        D(i:m,1)=real(temp);
        D(i:m,2)=imag(temp);
    end
    E(1:n)=D(1:n,1)+1i*D(1:n,2);
    H=triu(H,-1);
    for bella=1:n-1
        if(abs(H(bella+1,bella))<tol)
            H(bella+1,bella)=0;
        end
    end
end
