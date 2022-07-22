% HEAT EQUATION BACKWARD-DIFFERENCE 
%
clear all;   % b is right endpoint of time interval
m=10;n=50; L=1; b=.5; c=1; t=0; 
h=L/m; k=b/n;
r=c^2*k/(h^2);
 for i=1:m-1 
 w(i) = f2(i*h);
 end;
 w=w';
 for i=1:m-1
     a(i,i)=1+2*r;
     if i<m-1
        a(i,i+1)=-r;
        a(i+1,i)=-r;
     end
 end
 while t<b
    ww=a\w;
    w=ww;
    
    plot(w); hold on
    t=t+k;
 end
 w;