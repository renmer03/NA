clear all
a=0;b=1.6;N=16;h=(b-a)/N;x0=1;
w=zeros(N+1);w(1)=x0;
t=zeros(N+1);t(1)=a;
for i=1:N
    w(i+1)=w(i)+f(t(i),w(i))*h;
    t(i+1)=a+i*h; %t(i)+h
end

plot (t,w,'o'); hold on %axis equal;
y=t.^2+2*t+exp(t); plot(t,y)% exact solution
legend('euler solution','exact solution')
err=y-w;
