% WAVE EQUATION:  
%
 clear all;  % hold off
 L=1; T=1; c=4; m=10; n=20;
 h = L/m;  k = T/n; r = 0.5*k*c/h;
 for j = 2 : n+1
    w(1,j) = 0;
    w(m+1,j) = 0;
 end;
 w(1,1) = fw(0); w(m+1,1) = fw(L);
 for i = 2 : m 
    w(i,1) = fw((i-1)*h);
    w(i,2) = (1-r^2)*fw((i-1)*h)+0.5*r^2*(fw(i*h) ...
             +fw((i-2)*h))+k*g((i-1)*h);
 end;
 for j = 2 : n 
    for i = 2 : m 
       w(i,j+1) = 2*(1-r^2)*w(i,j)+r^2*(w(i+1,j) ...
                  +w(i-1,j))-w(i,j-1);
    end;
 end;
 fprintf('  I    X(I)     W(X(I),%2.1d)\n', T);
 
 for i = 1 : m+1 
    x(i) = (i-1)*h;
    fprintf('%3d %11.8f %13.8f\n', i, x(i), w(i,n+1));
 end;
%
%
axis([0 11 -1 1])
for i=1:n+1
  if i<=n/2
    plot(w(:,i));hold on; pause(0.1);
 else 
    plot(w(:,i),'r');pause(0.2);
 end
end
%
%
%  MOVIE
hold off
for k=1:n+1
  plot(fft(eye(k+16)))
    plot(w(:,k))
  axis ([0 11 -1 1])
  M(k)=getframe;
end
movie(M,10)
%
%EXACT SOLUTION
% xx=linspace(0,1);tt=0.5;
% u=sin(pi*xx).*cos(2*pi*tt);
