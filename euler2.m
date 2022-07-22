clear all
a=0;b=180;N=60;h=(b-a)/N;x0=293;
w=zeros(N+1);w(1)=x0;
t=zeros(N+1);t(1)=a;
for i=1:N
    w(i+1)=w(i)+f(t(i),w(i))*h+0.5*fp(t(i),w(i))*h^2;
    t(i+1)=a+i*h; %t(i)+h
end

plot (t,w,'o'); hold on %axis equal;
% axis([-0.01 1.6 0 11])
y=275+10*exp(-0.012*t); plot(t,y);% exact solution
% legend('order 2 -t=Taylor','Exact solution')
% err=y-w;
