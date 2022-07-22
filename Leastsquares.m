clear all
axis equal
x=[.78 .8 .82 .84 .86 .88]';
y=[.7 .72 .73 .74 .77 .77]';
m=length(x); plot (x,y,'*r'); hold on
p=sum(x); q=sum(y);r=x'*y;s=norm(x)^2; d=m*s-p^2;
a=(m*r-q*p)/d; b=(s*q-p*r)/d;
line=a*x+b;plot(x,line)