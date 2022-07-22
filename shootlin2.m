% Linear Shooting Algorithm 
%
 clear all;
 a=1; b=2; alpha=1; beta=2; n=20;
 fprintf('    t(i)         w(1,i)       w(2,i)\n');
 h = (b-a)/n; u1 = alpha; u2 = 0; v1 = 0; v2 = 1;
 for i = 1 : n 
     t = a+(i-1)*h;
     k11 = h*u2;
     k12 = h*(p(t)*u2+q(t)*u1+r(t));
     k21 = h*(u2+.5*k12);
     k22 = h*(p(t+.5*h)*(u2+0.5*k12)+q(t+.5*h)*(u1+0.5*k11)+r(t+.5*h));
     k31 = h*(u2+.5*k22);
     k32 = h*(p(t+.5*h)*(u2+.5*k22)+q(t+.5*h)*(u1+.5*k21)+rl(t+.5*h));
     k41 = h*(u2+k32);
     k42 = h*(p(t+h)*(u2+k32)+q(t+h)*(u1+k31)+r(t+h));
     u1 = u1+(k11+2*(k21+k31)+k41)/6;
     u2 = u2+(k12+2*(k22+k32)+k42)/6;
     %
     k11 = h*v2;
     k12 = h*(p(t)*v2+q(t)*v1);
     k21 = h*(v2+.5*k12);
     k22 = h*(p(t+.5*h)*(v2+.5*k12)+q(t+.5*h)*(v1+.5*k11));
     k31 = h*(v2+.5*k22);
     k32 = h*(p(t+.5*h)*(v2+.5*k22)+q(t+.5*h)*(v1+.5*k21));
     k41 = h*(v2+k32);
     k42 = h*(p(t+h)*(v2+k32)+q(t+h)*(v1+k31));
     v1 = v1+(k11+2*(k21+k31)+k41)/6;
     v2 = v2+(k12+2*(k22+k32)+k42)/6;
     u(1,i) = u1;   % u(1,:) is solution to first IVP
     u(2,i) = u2;
     v(1,i) = v1;   % v(1,:) is solution to second IVP 
     v(2,i) = v2;
 end;
 w1 = alpha;
 z = (beta-u(1,n))/v(1,n);
 t = a;
 i = 0;
 fprintf('%11.8f %11.8f %11.8f\n', t, w1, z);
 for i = 1 : n 
     t = a+i*h;
     w1 = u(1,i)+z*v(1,i);  % w1 IS SOLUTION TO BVP
     w2 = u(2,i)+z*v(2,i);  % w2 is derivative of w1
     fprintf('%11.8f %11.8f %11.8f\n', t, w1, w2);
 end;
