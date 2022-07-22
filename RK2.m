clear all
a=0;b=1.6;N=16;h=(b-a)/N;x0=1;
w=zeros(N+1);w(1)=x0;
t=zeros(N+1);t(1)=a;
for i=1:N
    k1=h*f(t(i),w(i));
    k2=h*f(t(i)+h,w(i)+k1);
    w(i+1)=w(i)+0.5*(k1+k2);
    t(i+1)=a+i*h; %t(i)+h
end

plot (t,w,'o'); hold on %axis equal;
axis([-0.01 1.6 0 11]);
y=t.^2+2*t+exp(t); plot(t,y);% exact solution
%legend('order 2 -t Taylor','Exact solution');
err=y-w;
