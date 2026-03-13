clear all
close all

sigma=1e-3;
b=5;
beta=0.1;
delta=1;

n=101;		% # spatial nodes
h=1/(n-1);
x=[0:h:1]';
uu=zeros(n,1);
A=2*diag(ones(n,1))-diag(ones(n-1,1),1)-diag(ones(n-1,1),-1);
A(1,1)=1;
A(n,n)=1;
A=A/h^2;
M=eye(n);
M(1,1)=1/2;
M(n,n)=1/2;

t=0;
dt=0.05;	% timestep
Tmax=50;

A=sigma*A+1/dt*M;

Iap=10;
bb1=zeros(n,1);
for i=1:n
   if(x(i)<=0.04)
     bb1(i)=Iap;
   end
end

k=1;
time(k)=t;
v_1(k)=uu(51);

while(t<=Tmax)
    t=t+dt;

    bb=zeros(n,1);
    if(t<=1) 
     bb=bb+bb1;
    end    
    bb=bb+b*uu.*(uu-beta).*(delta-uu);
    bb=M*(1/dt*uu+bb);

    uu_new=A\bb;

    uu=uu_new;
    
    k=k+1;
    time(k)=t;
    v_1(k)=uu(51);

    plot(x,uu,'or-')
    axis([0 1 -0.1 2])
    hold off 
    pause(0.001)
end

figure(2)
plot(time,v_1)
axis([0 Tmax -0.1 1.1])
