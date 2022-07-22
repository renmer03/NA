clear all;clf;
x=[-1 -075 -.5 .025 0 .25 .5 .75 1]';
y=[9, 7 5 4.5 5 6 8 11.5 15]';
m=length(x); plot(x,y,'*r'); hold on
n=2; T0=[0 0 1]; T1=[0 1 0];T2=[2 0-1];
A=zeros(n+1,n+1); b=zeros(n+1,1);
A(1,1)=n;A(1,2)=sum(x);A(1,3)=sum(polyval(T2,x));
A(2,1)=A(1,2); A(2,2)=x'*x;
A(2,3)=polyval(T2,x)'*x;
A(3,1)=A(1,3);A(3,2)=A(2,3);A(3,3)=polyval(T2,x)'*polyval(T2,x);
b(1)=sum(y);b(2)=x'*y;b(3)=y'*polyval(T2,x);
C=A\b; % Solve normal equations.
Curve=C(1)'*T0+C(2)'*T1+C(3)'*T2;
xx=linespace(-1,1); plot(xx,polyval(Curve,xx))