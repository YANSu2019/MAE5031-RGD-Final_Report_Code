%% Geometric Parameters
h0=1e-7;
alpha=0.005;
L=100;

%% Flow Parameters
mu=1.7e-5;
kB=1.3806505e-23;
m=40e-3/6.02e23;
U=25;
p0=1e5;
T=300;
Vm=sqrt(0.5*pi*kB*T/m);

%% main program

%% Calculate initialized Pressure
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

p2=p;

%% Iteration for Pressure
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
    
    Kn=zeros(1,xDim);
    dKn=zeros(1,xDim);
    for i=1:xDim
          Kn(i)=mu*Vm/p(i)/h(i)/h0/p0;
    end
    
    F=-0.9131.*Kn.*Kn.*Kn+2.8821.*Kn.*Kn+7.2833.*Kn+0.9644;
    
    dp=zeros(1,xDim);
    for i=2:xDim-1
        dp(i)=0.5*(p(i+1)-p(i-1))/dx;
    end
    
    dp(1)=(p(2)-p(1))/dx;
    dp(xDim)=(p(xDim)-p(xDim-1))/dx;
    
    for i=1:xDim
        dKn(i)=-(h(i)*dp(i)+p(i)*dh(i))*mu*Vm/p(i)/p(i)/h(i)/h(i)/h0/p0;
    end
    
    dF=(-3.*0.9131.*Kn.*Kn+2.*2.8821.*Kn+7.2833).*dKn;
    
    for i=2:xDim-1
          C(i,i-1)=-1.5*h(i)*h(i)*dh(i)*p(i)*dx*F(i)-0.5*h(i)*h(i)*h(i)*dF(i)*dx*p(i)+h(i)*h(i)*h(i)*p(i)*F(i)+3*mu*U*h(i)*dx/h0/p0;
          C(i,i)=-2*h(i)*h(i)*h(i)*p(i)*F(i)-6*mu*U*dh(i)*dx*dx/h0/p0;
          C(i,i+1)=1.5*h(i)*h(i)*dh(i)*p(i)*dx*F(i)+0.5*h(i)*h(i)*h(i)*dF(i)*dx*p(i)+h(i)*h(i)*h(i)*p(i)*F(i)-3*mu*U*h(i)*dx/h0/p0;
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


plot(p);