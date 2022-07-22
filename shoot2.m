% NONLINEAR SHOOTING ALGORITHM 
%
%      x'' = f(t,x,x'), a<=t<=b, x(a) = alpha, x(b) = beta:
%
 clear all  
 fprintf('  i        t(i)         W1(i)  \n');
 a=1; b=3; alpha=17; beta=43/3; m=8; n=20; tol = 1.e-5;
 h=(b-a)/n; s=(beta-alpha)/(b-a); 
 flag=0; k = 1;
 while k <= m && flag == 0
    w1(1) = alpha;
    w2(1) = s;
    u1 = 0 ;
    u2 = 1;
    for i = 1 : n 
       t = a+(i-1)*h;   t1 = t+0.5*h;
       k11 = h*w2(i);
       k12 = h*ft(t,w1(i),w2(i));
       k21 = h*(w2(i)+0.5*k12);
       k22 = h*ft(t1,w1(i)+0.5*k11,w2(i)+0.5*k12);
       k31 = h*(w2(i)+0.5*k22);
       k32 = h*ft(t1,w1(i)+0.5*k21,w2(i)+0.5*k22);
       k41 = h*(w2(i)+k32);
       k42 = h*ft(t+h,w1(i)+k31,w2(i)+k32);
       w1(i+1) = w1(i)+(k11+2*(k21+k31)+k41)/6;
       w2(i+1) = w2(i)+(k12+2*(k22+k32)+k42)/6;
       k11 = h*u2;
       k12 = h*(fy(t,w1(i),w2(i))*u1+fyp(t,w1(i),w2(i))*u2);
       k21 = h*(u2+0.5*k12);
       k22 = h*(fy(t1,w1(i),w2(i))*(u1+0.5*k11) ...
             +fyp(t1,w1(i),w2(i))*(u2+0.5*k21));
       k31 = h*(u2+0.5*k22);
       k32 = h*(fy(t1,w1(i),w2(i))*(u1+0.5*k21) ...
             +fyp(t1,w1(i),w2(i))*(u2+0.5*k22));
       k41 = h*(u2+k32);
       k42 = h*(fy(t+h,w1(i),w2(i))*(u1+k31) ...
             +fyp(t+h,w1(i),w2(i))*(u2+k32));
       u1  = u1+(k11+2*(k21+k31)+k41)/6;
       u2  = u2+(k12+2*(k22+k32)+k42)/6;
    end;
    if abs(w1(n+1)-beta) < tol 
       flag=1;
       i = 0;
       fprintf('%3d %13.8f %13.8f \n', i, a, alpha);
       for i = 1 : n 
          j = i+1;
          t = a+i*h;
          fprintf('%3d %13.8f %13.8f\n', i, t, w1(j));
       end;
%       fprintf(' t = %14.7e\n', tk);
    else
%      Improve shooting angle using Newton's
       s = s-(w1(n+1)-beta)/u1;
       k = k+1;
       if k > m 
	     flag=1; 
	       fprintf('iterations: %3d\n', k)
       end
    end;
 end;
