%The sextuple shift QR iteration
%input:a Hessenberg matrix
%output:after one step QR iteration
%W is 7*n matrix, ith column
function [W,H]=sextuple_shift_QR_iteration (H)
[n,~]=size(H);
a=poly(H(n-5:n,n-5:n));
temp=polyvalm(a,H(1:7,1:7));
vector=temp(1:7,1);
W=zeros(7,n-1);
for k=0:n-7
    w=household(vector);
    q=max([1,k]);
    r=min([k+8,n]);
    H(k+1:k+7,q:n)=H(k+1:k+7,q:n)-2*w*(w'*H(k+1:k+7,q:n));
    H(1:r,k+1:k+7)=H(1:r,k+1:k+7)-2*(H(1:r,k+1:k+7)*w)*w';
    W(1:7,k+1)=w;
    if(k<n-7)
        vector=H(k+2:k+8,k+1);
    end 
end
for i=1:5
    w=household(H(n-6+i:n,n-7+i));
    W(1:7-i,n-6+i)=w;
    H(n-6+i:n,n-7+i:n)=H(n-6+i:n,n-7+i:n)-2*w*(w'*H(n-6+i:n,n-7+i:n));
    H(1:n,n-6+i:n)=H(1:n,n-6+i:n)-2*(H(1:n,n-6+i:n)*w)*w';
end