mu=1.7e-5;
U=25;
h0=1e-7;
alpha=0.005;
L=100;
p0=1e5;
xDim=100;
x=linspace(0,L,xDim);
dx=x(2)-x(1);

h=zeros(1,xDim);
dh=zeros(1,xDim);
h=1-(x-100)*tan(alpha);
for i=1:xDim
    dh(i)=-tan(alpha);
end

A=zeros(xDim,xDim);
A(1,1)=1;
A(xDim,xDim)=1;
for i=2:xDim-1
    A(i,i-1)=h(i)*h(i)*h(i)-1.5*h(i)*h(i)*dh(i);
    A(i,i)=-2*h(i)*h(i)*h(i);
    A(i,i+1)=h(i)*h(i)*h(i)+1.5*h(i)*h(i)*dh(i);
end

p=ones(1,xDim);
cc=6*mu*U*dx*dx*dh(1)/(h0*p0);
B=cc*ones(1,xDim);
B(1)=1;
B(xDim)=1;

p=inv(A)*B';

%plot(p);
max(p)
%hold on

%%

p1=zeros(1,xDim);
C=zeros(xDim,xDim);
D=zeros(1,xDim);
C(1,1)=1;
C(xDim,xDim)=1;
err=zeros(1,xDim);
for i=1:xDim
    err(i)=abs(p1(i)-p(i));
end
ferr=max(err);

tstep=0;

while ferr>1e-7
    
    for i=2:xDim-1
          C(i,i-1)=-1.5*h(i)*h(i)*dh(i)*p(i)*dx+h(i)*h(i)*h(i)*p(i)+3*mu*U*h(i)*dx/h0/p0;
          C(i,i)=-2*h(i)*h(i)*h(i)*p(i)-6*mu*U*dh(i)*dx*dx/h0/p0;
          C(i,i+1)=1.5*h(i)*h(i)*dh(i)*p(i)*dx+h(i)*h(i)*h(i)*p(i)-3*mu*U*h(i)*dx/h0/p0;
          D(i)=-h(i)*h(i)*h(i)*0.25*(p(i+1)-p(i-1))*(p(i+1)-p(i-1));
    end
    C(1,1)=1;
    C(xDim,xDim)=1;
    D(1)=1;
    D(xDim)=1;
    p1=inv(C)*D';

   
    

    for i=1:xDim
         err(i)=abs(p1(i)-p(i));
    end
    ferr=max(err);
    p=p1;
    tstep=tstep+1;
end

p3=p;
plot(p);