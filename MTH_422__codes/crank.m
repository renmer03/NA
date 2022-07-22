% Crank-Nicolson Algorithm  
%
 clear all;
 m=10; L=1; T=0.5; n=50; c=1;
 h = L/m; k = T/n;
 r = c^2*k/(h^2);
 v(m) = 0;
 for i = 1 : m-1
    v(i) = f2(i*h);
 end;
%              Solving Tridiagonal Linear System
 l(1) = 1+r;
 u(1) = -r/(2*l(1));
 for i = 2 : m-2
    l(i) = 1+r+r*u(i-1)/2;
    u(i) = -r/(2*l(i));
 end;
 l(m-1) = 1+r+0.5*r*u(m-2);
 for j = 1 : n
    t(j) = j*k;
    z(1) = ((1-r)*v(1)+r*v(2)/2)/l(1);
    for i = 2 : m-1 
        z(i) = ((1-r)*v(i)+0.5*r*(v(i+1)+v(i-1)+z(i-1)))/l(i);
    end;
    v(m-1) = z(m-1);
    for i1 = 1 : m-2 
        i = m-i1-1;
        v(i) = z(i)-u(i)*v(i+1);
    end;
 end;
 fprintf('  i     x(i)       w(x(i),%5.2f)\n', T);
 fprintf('%3d %11.8f %13.8f\n', 0, 0, 0);
 %h
 for i = 1 : m-1 
    x(i) = i*h;
    fprintf('%3d %11.8f %13.8f\n', i, x(i), v(i));
 end;
 fprintf('%3d %11.8f %13.8f\n', 10, 1, 0);
x(m)=m*h;
