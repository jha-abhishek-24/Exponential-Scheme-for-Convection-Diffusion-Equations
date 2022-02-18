%Fourth order exponential scheme with TVD RK-3
%for solving 1-D convection diffusion IBVP
clc;
close all;
clear all;
%Space and time step sizes
h=0.04;
k=.05;
%Final time
time=20;
% Coefficients in the PDE
a=.01;
p=0.1;
l1=0;
l2=1;
n=((l2-l1)/h)+1;
x=linspace(l1,l2,n);
alpha=p*h/2*coth(p*h/(2*a));
alpha1=(a-alpha)/p;
alpha2=(a*(a-alpha))/p^2+h^2/6;
%Initial Condition
u00=zeros(1,n);
for i=1:n
    u00(1,i)=exp(5*x(i))*sin(pi*x(i));
end
%Boundary Conditions
g1=0;
g2=0;
% Numerical Scheme
A=zeros(n-2,n-2);
A(1,1)=1-(2*alpha2/h^2);
A(n-2,n-2)=1-(2*alpha2/h^2);
A(1,2)=alpha2/h^2+alpha1/(2*h);
A(n-2,n-3)=alpha2/h^2-alpha1/(2*h);
for i=2:n-3
    A(i,i)=1-(2*alpha2/h^2);
    A(i,i+1)=alpha2/h^2+alpha1/(2*h);
    A(i,i-1)=alpha2/h^2-alpha1/(2*h);
end
B=zeros(n-2,n-2);
B(1,1)=-2*alpha/h^2;
B(1,2)=alpha/h^2-p/(2*h);
B(n-2,n-3)=alpha/h^2+p/(2*h);
B(n-2,n-2)=-2*alpha/h^2;
for i=2:n-3
    B(i,i)=-2*alpha/h^2;
    B(i,i+1)=alpha/h^2-p/(2*h);
    B(i,i-1)=alpha/h^2+p/(2*h);
end
    g=zeros(n-2,1);
g(1,1)=(alpha/h^2+p/(2*h))*g1;
g(n-2,1)=(alpha/h^2-p/(2*h))*g2;
un=zeros(n-2,1);
u0=zeros(n-2,1);
u0=u00(2:n-1)';
t=0;
 while t<=time
 L1=(A\B)*(u0)+(A\g);
 C1=u0+k*L1;
 L2=(A\B)*(C1)+(A\g);
 C2=(3/4)*u0+(1/4)*C1+(1/4)*k*L2;
 L3=(A\B)*(C2)+(A\g);
 un=(1/3)*u0+(2/3)*C2+(2/3)*k*L3;
 u0=un;
 t=t+k;
 end
una=zeros(n,1);
una(1)=g1;
una(n)=g2;
una(2:n-1)=un(1:n-2);
%Analytical Solution
ua=zeros(n,1);
for i=1:n
    ua(i,1)=exp(5*x(i)-time*(0.01*pi^2+0.25))*sin(pi*x(i));
end
%Error
err=norm((ua-una),2)
%Order of convergence
order=log2(7.6825e-05/6.3298e-06) 
plot(x,una,'*b-');
hold on
plot(x,ua,'*r');
axis([0 1 0 0.035])
xlabel('x');
ylabel('u(x,20)');
legend('una', 'ua');
title('Solution with h=0.04 and time, T=20')
%------------------------------------------------------------




