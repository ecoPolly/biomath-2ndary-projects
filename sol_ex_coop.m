% time interval
t0=0;
tmax=10;

% discretization parameter
h=0.001;

% enter paraters
s0=5;
e0=1;
k2=1;
k4=2;

coop=input('coop = (positive 1, none 0, negative -1) ');
switch(coop)
   case 1
      disp('positive cooperativity')
      km1=500;
      k1=0.501;
      km3=500;
      k3=502000;
   case 0
      disp('no cooperativity')
      km1=50;
      k1=102;
      km3=50;
      k3=26;
   case -1
      disp('negative cooperativity')
      km1=500;
      k1=1002;
      km3=500;
      k3=5.02;
end

param(1)=e0;
param(2)=k2;
param(3)=k4;
param(4)=k1;
param(5)=km1;
param(6)=k3;
param(7)=km3;

% solution
options=odeset('RelTol',5.e-13 ,'AbsTol',1.e-13*ones(1,3),'InitialStep',1.e-5);
[T,Y]=ode15s(@eq_fun,[t0:h:tmax],[s0 0 0],options,param);

options=odeset('RelTol',5.e-13 ,'AbsTol',1.e-13,'InitialStep',1.e-5);
[T,Y_qs]=ode15s(@eq_u_fun,[t0:h:tmax],s0,options,param);

s_ex=Y(:,1);
c1_ex=Y(:,2);
c2_ex=Y(:,3);

s_qs=Y_qs;

K_c=(km3+k4)/k3;
K_m=(km1+k2)/k1;

c1_qs=K_c*e0*s_qs./(K_c*s_qs+s_qs.^2+K_m*K_c);
c2_qs=e0*s_qs.^2./(K_c*s_qs+s_qs.^2+K_m*K_c);

% compute relative error between Y and Y_qs
err_s=norm(s_qs-s_ex,inf)/norm(s_ex,inf);
err_c1=norm(c1_qs-c1_ex,inf)/norm(c1_ex,inf);
err_c2=norm(c2_qs-c2_ex,inf)/norm(c2_ex,inf);
disp(['error on component s = ',num2str(err_s)]);
disp(['error on component c1 = ',num2str(err_c1)]);
disp(['error on component c2 = ',num2str(err_c2)]);

% plot
V_ex=k2*c1_ex+k4*c2_ex;
V_qs=k2*c1_qs+k4*c2_qs;

plot(s_ex,V_ex,'b','Linewidth',2)
hold on
plot(s_qs,V_qs,'r','Linewidth',2)
title('V_{ex} (blue), V_{qs} (red)','Fontsize',24)


% function: eq_fun, exact dynamics
function dy=eq_fun(t,y,param)

e0=param(1);
k2=param(2);
k4=param(3);
k1=param(4);
km1=param(5);
k3=param(6);
km3=param(7);

s=y(1);
c1=y(2);
c2=y(3);

dy=zeros(3,1);

e=e0-c1-c2;

dy(1)=-k1*s*e+km1*c1-k3*s*c1+km3*c2;
dy(2)=k1*s*e-(km1+k2)*c1-k3*s*c1+(k4+km3)*c2;
dy(3)=k3*s*c1-(k4+km3)*c2;

end


% function: eq_qs_fun, uniform approximation
function dy=eq_u_fun(t,y,param)

e0=param(1);
k2=param(2);
k4=param(3);
k1=param(4);
km1=param(5);
k3=param(6);
km3=param(7);

s=y(1);

K_c=(km3+k4)/k3;
K_m=(km1+k2)/k1;

c1=K_c*e0*s./(K_c*s+s.^2+K_m*K_c);
c2=e0*s.^2./(K_c*s+s.^2+K_m*K_c);

dy=-(k2*c1+k4*c2);

end
