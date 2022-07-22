clear all
a=0;b=1.6;N=16;h=(b-a)/N;x0=1;
w=zeros(N+1);w(1)=x0;
t=zeros(N+1);t(1)=a;
for i=1:N
    w(i+1)=w(i)+f(t(i),w(i))*h+0.5*fp(t(i),w(i))*h^2+1/6*fp2(t(i),w(i))*h^3+1/24*fp3(t(i),w(i))*h^4;
    t(i+1)=a+i*h; %t(i)+h
end

plot (t,w,'o','MarkerFaceColor','r'); hold on %axis equal;
axis([-0.01 1.6 0 11])
y=t.^2+2*t+exp(t); plot(t,y);% exact solution
legend('order 4 -t=Taylor','Exact solution')
err=y-w;
