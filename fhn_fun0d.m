function dydt=rm_fun(t,y,Iap,Tsti)

v=y(1);
w=y(2);

% PARAMETERS
G=1.5;
eta1=4.4;
eta2=0.012;
vth=13;
vp=100;

if(t<=Tsti)
  Iapp=Iap;
else
  Iapp=0;
end

dydt(1)=-G*v.*(1-v/vth).*(1-v/vp)-eta1*v.*w+Iapp;
dydt(2)=eta2*(v/vp-w);
dydt=dydt';
