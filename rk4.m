% Numerical Analysis
% by Jorge Rebaza
% RUNGE-KUTTA (ORDER 4) 
% INPUT:   ENDPOINTS A,B; INITIAL CONDITION X0; INTEGER N.
%
clear all; %hold off;
a=0; b=1.6; N=16; h = (b-a)/N; x0=[1,1,-1]; m = 3;
w=x0;
t=a;
for i = 1:N 
    k1 = h*f(t,w(1),w(2),w(3));
    k2 = h*f(t+h/2, w(1)+k1(1)/2,w(2)+ k1(2)/2, w(3)+k1(3)/2);
    k3 = h*f(t+h/2, w(1)+k2(1)/2,w(2)+ k2(2)/2, w(3)+k2(3)/2);
    k4 = h*f(t+h,w(1)+k3(1)/2,w(2)+ k3(2)/2, w(3)+k3(3)/2);
    w = w +(k1+2*k2+2*k3+k4)/6;
    t= a+i*h; % = t(i)+h
    fprintf('%5.3f', t);
    for j=1:m
        fprintf('%11.8f', w(j));
    end
    fprintf('\n');
end;
% plot(t,w,'o','MarkerFaceColor','r'); hold on
% axis([-0.01 1.6 0 11])
% y=t.^2+2*t+exp(t); plot(t,y);% Exact solution
% legend('Runge-Kutta order 4','Exact Solution')
% err=y-w;
% [y w abs(y-w)]