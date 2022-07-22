% % Numerical Methods
%
% PARABOLIC EQ. FORWARD-DIFFERENCE 
%
clear all;   % b is right endpoint of time interval
m=16;n=125; L=1; b=0.125; c=1; t=0;   % n=125 stable, n=40 unstable
%
h=L/m; k=0.002; %k=b/n;
r=c^2*k/(h^2);
a=zeros(m-1,m-1);
for i=1:m-1 
   w(i) = f3(i*h);
end;
w=w';
%
% Form matrix a(i,j) here
%
 while t<b   %% 
    ww=a*w;
    w=ww; v=[0 w' 0]; 
    %plot(v); hold on
    t=t+k;
 end
 % comparing with exact solution
 tt=linspace(0,1,m+1);
 plot(tt,v,'ro'); hold on
%  y=exp(-1).*sin(pi*tt)+exp(-9)*sin(3*pi*tt);plot(tt,y) % EXACT at t=1
 