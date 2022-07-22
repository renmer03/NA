% Elliptic  Equation
%
 clear all;
 a=0; b=2; c=0; d=1; n=10; m=10; max=200; flag=0; TOL=1.e-6;
 h = (b-a)/n;  k = (d-c)/m;
 for i = 0 : n
    x(i+1) = a+i*h;
 end;
 for j = 0 : m
    y(j+1) = c+j*k;
 end;
 for i = 1 : n-1
    w(i+1,1) = gp(x(i+1),y(1));
    w(i+1,m+1) = gp(x(i+1),y(m+1));
 end;
 for j = 0 : m
    w(1,j+1) = gp(x(1),y(j+1));
    w(n+1,j+1) = gp(x(n+1),y(j+1));
 end;
 lam = h*h/(k*k); mu = 2*(1+lam);
%
 L = 1;    % Solution of Linear System (Gauss-Seidel) until line 85
 while L <= max && flag == 0
     z = (-h*h*fp(x(2),y(m)) + gp(a,y(m)) + lam*gp(x(2),d) + lam*w(2,m-1) + w(3,m))/mu;
     e = abs( w(2,m)-z); w(2,m) = z;
     for i = 2 : n-2
        z = (-h*h*fp(x(i+1),y(m)) + lam*gp(x(i+1),d)+w(i,m) + w(i+2,m) + lam*w(i+1,m-1))/mu;
        if abs(w(i+1,m)-z) > e
           e = abs( w(i+1,m) - z );
        end;
        w(i+1,m) = z;
     end;
     z = (-h*h*fp(x(n),y(m)) + gp(b,y(m)) + lam*gp(x(n),d) + w(n-1,m) + lam*w(n,m-1))/mu;
     if abs( w(n,m)-z) > e
        e = abs( w(n,m)-z);
     end;
     w(n,m) = z;
     for k = 2 : m-2
        j = m-k;
        z = (-h*h*fp(x(2),y(j+1)) + gp(a,y(j+1)) + lam*w(2,j+2) + lam*w(2,j) + w(3,j+1))/mu;
        if abs(w(2,j+1)-z) > e
           e = abs(w(2,j+1)-z);
        end;
        w(2,j+1) = z;
        for i = 2 : n-2
           z = (-h*h*fp(x(i+1),y(j+1)) + w(i,j+1) + lam*w(i+1,j+2) + lam*w(i+1,j) + w(i+2,j+1))/mu;
           if abs(w(i+1,j+1)-z) > e
              e = abs(w(i+1,j+1)-z);
           end;
           w(i+1,j+1) = z;
        end;
        z = (-h*h*fp(x(n),y(j+1))+gp(b,y(j+1))+w(n-1,j+1)+lam*w(n,j+2)+lam*w(n,j))/mu;
        if abs(w(n,j+1)-z) > e
           e = abs(w(n,j+1)-z);
        end;
        w(n,j+1) = z;
     end;
     z = (-h*h*fp(x(2),y(2)) + lam*gp(x(2),c) + gp(a,y(2)) + lam*w(2,3) + w(3,2))/mu;
     if abs(w(2,2)-z) > e
        e = abs(w(2,2)-z);
     end;
     w(2,2) = z;
     for i = 2 : n-2
        z = (-h*h*fp(x(i+1),y(2)) + lam*gp(x(i+1),c) + w(i+2,2) + w(i,2) + lam*w(i+1,3))/mu;
        if abs(w(i+1,2)-z) > e
           e = abs(w(i+1,2)-z);
        end;
        w(i+1,2) = z;
     end;
     z = (-h*h*fp(x(n),y(2)) + lam*gp(x(n),c) + gp(b,y(2)) + w(n-1,2) + lam*w(n,3))/mu;
     if abs(w(n,2)-z) > e
        e = abs(w(n,2)-z);
     end;
     w(n,2) = z;
     if e <= TOL
        flag=1;
        fprintf('  i  j    x(i)        y(j)         w(i,j)\n\n');
        for i = 1 : n
           for j = 1 : m 
              fprintf('%3d %2d %11.8f %11.8f %13.8f\n',i,j,x(i+1),y(j+1),w(i+1,j+1));
           end;
        end;
     else
        L = L+1;   % End of Gauss-Seidel
     end;
 end;
 if L > max 
    fprintf('Method fails after iteration number %d\n', max)
 end;
 surf(y,x,w)
