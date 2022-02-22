%Do household transformation (I-2ww')x=e
%input:vector x
%output:Household vector w and matrix H
%function [H,w,beta]=household(x)
%n=size(x,1);
%yita=norm(x,inf);
%x=x/yita;
%w=zeros(n,1);
%a=x(2:n,1)'*x(2:n,1);
%w(2:n,1)=x(2:n);
%if a==0
 %   beta=0;
%else
%    alpha=sqrt(w(1,1)^2+a);
%    if x(1)<0
%        w(1)=x(1)-alpha;
%    else
%        w(1)=-a/(x(1)+alpha);
%    end
%    beta=2*w(1,1)^2/(a+w(1,1)^2);
%    w=w/w(1,1);
%end
%H=eye(n);
%H=H-beta*(w*w');
%end
function [w]=household(x) 
    n=size(x,1); w=zeros(n,1); 
    a=norm(x(2:n),2); 
    b=norm(x,2); 
    if a==0
        w=0;
    else
        if x(1)<0 
            w(1)=x(1)-b; 
        else
            w(1)=-a^2/(x(1)+b);
        end
        for i=2:n 
            w(i)=x(i); 
        end
        t=norm(w,2);
        w=w/t; 
    end
end