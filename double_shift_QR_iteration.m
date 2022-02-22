%double shift QR iteration
%input:a Hessenberg matrix
%output:after one step QR iteration
%W is 3*n matrix, ith column
function [W,H]=double_shift_QR_iteration(H)
    [n,~]=size(H);
    s=H(n-1,n-1)+H(n,n);
    t=H(n-1,n-1)*H(n,n)-H(n-1,n)*H(n,n-1);
    x=H(1,1)*H(1,1)+H(1,2)*H(2,1)-s*H(1,1)+t;
    y=H(2,1)*(H(1,1)+H(2,2)-s);
    z=H(2,1)*H(3,2);
    vector=[x;y;z];
    W=zeros(3,n-1);
    for k=0:n-3
        w=household(vector);
        W(1:3,k+1)=w;
        q=max([1,k]);
        r=min([k+4,n]);
        H(k+1:k+3,q:n)=H(k+1:k+3,q:n)-2*w*(w'*H(k+1:k+3,q:n));
        H(1:r,k+1:k+3)=H(1:r,k+1:k+3)-2*(H(1:r,k+1:k+3)*w)*w';
        %Q(1:al,k+1:k+3)=Q(1:al,k+1:k+3)-2*(Q(1:al,k+1:k+3)*w)*w';
        x=H(k+2,k+1);
        y=H(k+3,k+1);
        if(k<n-3)
            z=H(k+4,k+1);
        end
        vector=[x;y;z];
    end
    w=household([x;y]);
    W(1:2,n-1)=w;
    H(n-1:n,n-2:n)=H(n-1:n,n-2:n)-2*w*(w'*H(n-1:n,n-2:n));
    H(1:n,n-1:n)=H(1:n,n-1:n)-2*(H(1:n,n-1:n)*w)*w';
end