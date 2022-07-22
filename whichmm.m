a=[1];
b=[1];
c=[1];
d=[1];
e=[1];
flag=0;

while flag==0
    n=size(a)
    for i=1:n
    a=a+1;
    F(i)=1/a(i)+1/b(i)+1/c(i)+1/d(i)+1/e(i)
    if F(i)==1
        F append F(i)
    else
        b(i)=b(i)+1;