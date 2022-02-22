%double_shift_QR_algorithm
%input:A dense matrix H
%output:A n*2 matrix with all of the eigenvalue of A;
%The first column is the real part and the second column is the second
%column is the imaginary part, which the a+bi is front of a-bi and sort by
%the length of eigenvalues.
function [E,H,Q]=sextuple_shift_QR_algorithm(A,flag)
    [n,~]=size(A);
    D=zeros(n,2);
    E=zeros(n,1);
    H=A;
    Q=eye(n);
    if flag==1
        %[Q,H]=hessenberg(H);
        [Q,H]=hess(H);
    end
    i=1;
    m=n;
    tol=1e-15;
    while(m-i+1>7)
        [W,H(i:m,i:m)]=sextuple_shift_QR_iteration(H(i:m,i:m));
        [~,avavaava]=size(W);
        for tp=1:avavaava-5
            w=W(1:7,tp);
            Q(1:n,i-1+tp:i+5+tp)=Q(1:n,i-1+tp:i+5+tp)-2*(Q(1:n,i-1+tp:i+5+tp)*w)*w';
            H(1:i-1,i-1+tp:i+5+tp)=H(1:i-1,i-1+tp:i+5+tp)-2*(H(1:i-1,i-1+tp:i+5+tp)*w)*w';
            H(i-1+tp:i+5+tp,m+1:n)=H(i-1+tp:i+5+tp,m+1:n)-2*w*(w'*H(i-1+tp:i+5+tp,m+1:n));
        end
        for ilp=1:5
            w=W(1:7-ilp,avavaava-5+ilp);
            Q(1:n,i-1+avavaava+ilp-5:i+avavaava)=Q(1:n,i-1+avavaava+ilp-5:i+avavaava)-2*(Q(1:n,i-1+avavaava+ilp-5:i+avavaava)*w)*w';
            H(1:i-1,i-1+avavaava+ilp-5:i+avavaava)=H(1:i-1,i-1+avavaava+ilp-5:i+avavaava)-2*(H(1:i-1,i-1+avavaava+ilp-5:i+avavaava)*w)*w';
            H(i-1+avavaava+ilp-5:i+avavaava,m+1:n)=H(i-1+avavaava+ilp-5:i+avavaava,m+1:n)-2*w*(w'*H(i-1+avavaava+ilp-5:i+avavaava,m+1:n));
        end
        while(m-i+1>7)
            py=0;
            pj=0;
            if(abs(H(m,m-1))<tol)
                py=1;
            end
            if(abs(H(m-1,m-2))<tol)
                py=2;
            end
            if(abs(H(m-2,m-3))<tol)
                py=3;
            end
            if(abs(H(m-3,m-4))<tol)
                py=4;
            end
            if(abs(H(m-4,m-5))<tol)
                py=5;
            end
            if(abs(H(m-5,m-6))<tol)
                py=6;
            end
            if(abs(H(i+1,i))<tol)
                pj=1;
            end
            if(abs(H(i+2,i+1))<tol)
                pj=2;
            end
            if(abs(H(i+3,i+2))<tol)
                pj=3;
            end
            if(abs(H(i+4,i+3))<tol)
                pj=4;
            end
            if(abs(H(i+5,i+4))<tol)
                pj=5;
            end
            if(abs(H(i+6,i+5))<tol)
                pj=6;
            end
            tlpo=0;
            if(py~=0)
                tlpo=tlpo+1;
            end
            if(pj~=0)
                tlpo=tlpo+1;
            end
            if(tlpo==0)
                break;
            end
            if(py~=0)
                switch py
                    case 1
                        H(m,m-1)=0;
                        D(m,1)=H(m,m);
                        m=m-1;
                    case 2
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
                        m=m-2;
                    case 3
                        [E0,H(m-2:m,m-2:m),Q0]=double_shift_QR_algorithm(H(m-2:m,m-2:m),0);
                        Q(1:n,m-2:m)=Q(1:n,m-2:m)*Q0;
                        H(m-2:m,m+1:n)=Q0'*H(m-2:m,m+1:n);
                        H(1:m-3,m-2:m)=H(1:m-3,m-2:m)*Q0;
                        D(m-2:m,1)=real(E0);
                        D(m-2:m,2)=imag(E0);
                        m=m-3;
                    case 4
                        [E0,H(m-3:m,m-3:m),Q0]=double_shift_QR_algorithm(H(m-3:m,m-3:m),0);
                        Q(1:n,m-3:m)=Q(1:n,m-3:m)*Q0;
                        H(m-3:m,m+1:n)=Q0'*H(m-3:m,m+1:n);
                        H(1:m-4,m-3:m)=H(1:m-4,m-3:m)*Q0;
                        D(m-3:m,1)=real(E0);
                        D(m-3:m,2)=imag(E0);
                        m=m-4;
                    case 5
                        [E0,H(m-4:m,m-4:m),Q0]=double_shift_QR_algorithm(H(m-4:m,m-4:m),0);
                        Q(1:n,m-4:m)=Q(1:n,m-4:m)*Q0;
                        H(m-4:m,m+1:n)=Q0'*H(m-4:m,m+1:n);
                        H(1:m-5,m-4:m)=H(1:m-5,m-4:m)*Q0;
                        D(m-4:m,1)=real(E0);
                        D(m-4:m,2)=imag(E0);
                        m=m-5;
                    case 6
                        [E0,H(m-5:m,m-5:m),Q0]=double_shift_QR_algorithm(H(m-5:m,m-5:m),0);
                        Q(1:n,m-5:m)=Q(1:n,m-5:m)*Q0;
                        H(m-5:m,m+1:n)=Q0'*H(m-5:m,m+1:n);
                        H(1:m-6,m-5:m)=H(1:m-6,m-5:m)*Q0;
                        D(m-5:m,1)=real(E0);
                        D(m-5:m,2)=imag(E0);
                        m=m-6;
                end
            end
            if (pj~=0)
                switch pj
                    case 1
                        H(i+1,i)=0;
                        D(i,1)=H(i,i);
                        i=i+1;
                    case 2
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
                        i=i+2;
                    case 3
                        [E0,H(i:i+2,i:i+2),Q0]=double_shift_QR_algorithm(H(i:i+2,i:i+2),0);
                        Q(1:n,i:i+2)=Q(1:n,i:i+2)*Q0;
                        H(i:i+2,i+3:n)=Q0'*H(i:i+2,i+3:n);
                        H(1:i-1,i:i+2)=H(1:i-1,i:i+2)*Q0;
                        D(i:i+2,1)=real(E0);
                        D(i:i+2,2)=imag(E0);
                        i=i+3;
                    case 4
                        [E0,H(i:i+3,i:i+3),Q0]=double_shift_QR_algorithm(H(i:i+3,i:i+3),0);
                        Q(1:n,i:i+3)=Q(1:n,i:i+3)*Q0;
                        H(i:i+3,i+4:n)=Q0'*H(i:i+3,i+4:n);
                        H(1:i-1,i:i+3)=H(1:i-1,i:i+3)*Q0;
                        D(i:i+3,1)=real(E0);
                        D(i:i+3,2)=imag(E0);
                        i=i+4;
                    case 5
                        [E0,H(i:i+4,i:i+4),Q0]=double_shift_QR_algorithm(H(i:i+4,i:i+4),0);
                        Q(1:n,i:i+4)=Q(1:n,i:i+4)*Q0;
                        H(i:i+4,i+5:n)=Q0'*H(i:i+4,i+5:n);
                        H(1:i-1,i:i+4)=H(1:i-1,i:i+4)*Q0;
                        D(i:i+4,1)=real(E0);
                        D(i:i+4,2)=imag(E0);
                        i=i+5;
                    case 6
                        [E0,H(i:i+5,i:i+5),Q0]=double_shift_QR_algorithm(H(i:i+5,i:i+5),0);
                        Q(1:n,i:i+5)=Q(1:n,i:i+5)*Q0;
                        H(i:i+5,i+6:n)=Q0'*H(i:i+5,i+6:n);
                        H(1:i-1,i:i+5)=H(1:i-1,i:i+5)*Q0;
                        D(i:i+5,1)=real(E0);
                        D(i:i+5,2)=imag(E0);
                        i=i+6;
                end
            end
        end
    end
    if(m>i)
        [temp,H(i:m,i:m),Q0]=double_shift_QR_algorithm(H(i:m,i:m),0);
        Q(1:n,i:m)=Q(1:n,i:m)*Q0;
        H(i:m,m+1:n)=Q0'*H(i:m,m+1:n);
        H(1:i-1,i:m)=H(1:i-1,i:m)*Q0;
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
