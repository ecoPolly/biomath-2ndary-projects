function [w1,w2,w3,bb,ina,ik]=HH_fun(w1,w2,w3,uu,dt,nno)

bb=zeros(nno,1);
ina=zeros(nno,1);
ik=zeros(nno,1);

% CONSTANTS
cm=1;
gna=120;
gk=36;
gl=0.3;
Vna=115;
Vk=-12;
Vl=10.6;

for i=1:nno
  V=uu(i);
  m=w1(i);
  h=w2(i);
  n=w3(i);

  % UPDATE GATING 
  am=0.1*(25-V)./(exp((25-V)/10)-1);
  bm=4*exp(-V/18);
  ah=0.07*exp(-V/20);
  bh=1./(exp((30-V)/10)+1);
  an=0.01*(10-V)./(exp((10-V)/10)-1);
  bn=0.125*exp(-V/80);

  m=(m+dt*am)/(1+dt*(am+bm));
  h=(h+dt*ah)/(1+dt*(ah+bh));
  n=(n+dt*an)/(1+dt*(an+bn));
  
  w1(i)=m;
  w2(i)=h;
  w3(i)=n;

  ina(i) = gna*m.^3.*h.*(V-Vna);
  ik(i) = gk*n.^4.*(V-Vk);
  ion=ina(i)+ik(i)+gl*(V-Vl);
  bb(i)=-ion;
  
end
