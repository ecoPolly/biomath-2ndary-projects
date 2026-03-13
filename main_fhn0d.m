% INITIAL CONDITIONS
v0=0;
w0=0;

y0=[v0 w0]';

% PARAMETERS
iap=input('iap = ');
Tsti=input('Tsti = ');

options=odeset('RelTol',5.e-13 ,'AbsTol',[1.e-13 1.e-13],'InitialStep',1.e-5,'MaxStep',5);

[t1,y1]=ode15s(@fhn_fun0d,[0 200],y0,options,iap,Tsti);

plot(t1,y1(:,1))


function dydt=fhn_fun0d(t,y,Iap,Tsti)

v=y(1);
w=y(2);

% PARAMETERS
b=5;
c=1;
beta=0.1;
delta=1;
gamma=0.25;
e=0.1;

if(t<=Tsti)
  Iapp=Iap;
else
  Iapp=0;
end

dydt(1)=b*v.*(v-beta).*(delta-v)-c*w+Iapp;
dydt(2)=e*(v-gamma*w);
dydt=dydt';

end
