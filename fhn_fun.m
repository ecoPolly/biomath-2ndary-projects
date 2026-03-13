function [ww,bb]=fhn_fun(ww,uu,dt,nno)

bb=zeros(nno,1);

% PARAMETERS
b=5;
c=1;
beta=0.1;
delta=1;
gamma=0.25;
e=0.1;

% LOOP ON NODES
for i=1:nno
	u=uu(i);
	w=ww(i);

	w=(w+dt*e*u)/(1+dt*e*gamma);
	Ion=b*u*(u-beta)*(delta-u)-c*w;

        ww(i)=w;
	bb(i)=Ion;
end	
