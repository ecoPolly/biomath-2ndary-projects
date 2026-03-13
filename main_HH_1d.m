close all;
clear all;

% INPUT
nno=101;
dt=0.05;
Tend=50;
h=1/(nno-1);
x=linspace(0,1,nno);

% INITIALIZE
y0=[2.7570e-04   5.2934e-02   5.9611e-01   3.1768e-01]';
uu=y0(1)*ones(nno,1);
ww1=y0(2)*ones(nno,1);
ww2=y0(3)*ones(nno,1);
ww3=y0(4)*ones(nno,1);

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
m_1(k)=ww1(51);
h_1(k)=ww2(51);
n_1(k)=ww3(51);

bb1=zeros(nno,1);
bb1(1:5)=Iap;

while(t<=Tend)
	t=t+dt;
	
        [ww1,ww2,ww3,bb,ina,ik]=HH_fun(ww1,ww2,ww3,uu,dt,nno);

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
        m_1(k)=ww1(51);
        h_1(k)=ww2(51);
        n_1(k)=ww3(51);

	% PLOT
	plot(x,uu,'o-r')
        axis([0 1 -20 120])
	pause(0.001)
	hold off
end


