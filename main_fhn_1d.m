close all;
clear all;

% INPUT
nno=101;
dt=0.02;
Tend=100;
h=1/(nno-1);
x=linspace(0,1,nno);

% INITIALIZE
y0=[0 0]';
uu=y0(1)*ones(nno,1);
ww=y0(2)*ones(nno,1);

% ASSEMBLE MATRIX
sigma=0.001;
one=ones(nno,1);
AA=spdiags([-one 2*one -one],-1:1,nno,nno);
AA(1,1)=1;
AA(nno,nno)=1;
one(1)=1/2;
one(nno)=1/2;
AA=spdiags(one,0,nno,nno)/dt+sigma/h^2*AA;

t=0;
k=1;
Iap = 100;
Tsti = 1;

time(k)=t;
v_1(k)=uu(51);
m_1(k)=ww(51);

bb1=zeros(nno,1);
bb1(1:5)=Iap;

while(t<=Tend)
	t=t+dt;
	
        [ww,bb]=fhn_fun(ww,uu,dt,nno);

	if(t<=Tsti)
	  bb=bb+bb1;
        end

	bb=bb+uu/dt;

	bb(1)=bb(1)/2;
	bb(nno)=bb(nno)/2;

	uu_new=AA\bb;
    
	% UPDATE
	uu=uu_new;

	k=k+1;
	time(k)=t;
        v_1(k)=uu(51);
        m_1(k)=ww(51);

	% PLOT
	plot(x,uu,'o-r')
        axis([0 1 -0.5 1.2])
	pause(0.001)
	hold off
end


