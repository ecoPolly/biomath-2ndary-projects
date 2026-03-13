function dydt=HH(t,y,di,Tsti,Iapp)

% VARIABLES
V=y(1);
m=y(2);
h=y(3);
n=y(4);

% CONSTANTS
cm=1;
gna=120;
gk=36;
gl=0.3;
Vna=115;
Vk=-12;
Vl=10.6;

% COMPUTE RHS
am=0.1*(25-V)./(exp((25-V)/10)-1);
bm=4*exp(-V/18);
ah=0.07*exp(-V/20);
bh=1./(exp((30-V)/10)+1);
an=0.01*(10-V)./(exp((10-V)/10)-1);
bn=0.125*exp(-V/80);

dydt(1)=-1/cm*(gna*m.^3.*h.*(V-Vna)+gk*n.^4.*(V-Vk)+gl*(V-Vl));
if(mod(t,di)<=Tsti)
  dydt(1)=dydt(1)+1/cm*Iapp;
end
dydt(2)=am.*(1-m)-bm.*m;
dydt(3)=ah.*(1-h)-bh.*h;
dydt(4)=an.*(1-n)-bn.*n;

dydt=dydt';
