% time interval
t0=0;
tmax=10;

% discretization parameter
h=0.001;

% enter paraters
s0=5;
i0=2;
e0=1;
km1=50;
k1=102;
k2=1;
km3=50;
k3=26;

param(1)=e0;
param(2)=k2;
param(3)=k1;
param(4)=km1;
param(5)=k3;
param(6)=km3;

% solution
options=odeset('RelTol',5.e-13 ,'AbsTol',1.e-13*ones(1,4),'InitialStep',1.e-5);
[T,Y]=ode15s(@eq_fun,[t0:h:tmax],[s0 i0 0 0],options,param);

options=odeset('RelTol',5.e-13 ,'AbsTol',1.e-13*ones(1,2),'InitialStep',1.e-5);
[T,Y_qs]=ode15s(@eq_qs_fun,[t0:h:tmax],[s0 i0],options,param);

s_ex=Y(:,1);
i_ex=Y(:,2);
c1_ex=Y(:,3);
c2_ex=Y(:,4);

s_qs=Y_qs(:,1);
i_qs=Y_qs(:,2);

K_i=km3/k3;
K_m=(km1+k2)/k1;

c1_qs=e0*s_qs./(s_qs+K_m*(1+i_qs/K_i));
c2_qs=K_m/K_i*e0*i_qs./(s_qs+K_m*(1+i_qs/K_i));

% compute relative error between Y and Y_qs
err_s=norm(s_qs-s_ex,inf)/norm(s_ex,inf);
err_i=norm(i_qs-i_ex,inf)/norm(i_ex,inf);
err_c1=norm(c1_qs-c1_ex,inf)/norm(c1_ex,inf);
err_c2=norm(c2_qs-c2_ex,inf)/norm(c2_ex,inf);
disp(['error on component s = ',num2str(err_s)]);
disp(['error on component i = ',num2str(err_i)]);
disp(['error on component c1 = ',num2str(err_c1)]);
disp(['error on component c2 = ',num2str(err_c2)]);

% plot
subplot(2,2,1)
plot(T,s_ex,'b','Linewidth',2)
hold on
plot(T,s_qs,'r','Linewidth',2)
title('s_{ex} (blue), s_{qs} (red)','Fontsize',24)

subplot(2,2,2)
plot(T,i_ex,'b','Linewidth',2)
hold on
plot(T,i_qs,'r','Linewidth',2)
title('i_{ex} (blue), i_{u} (red)','Fontsize',24)

subplot(2,2,3)
plot(T,c1_ex,'b','Linewidth',2)
hold on
plot(T,c1_qs,'r','Linewidth',2)
title('c1_{ex} (blue), c1_{qs} (red)','Fontsize',24)

subplot(2,2,4)
plot(T,c2_ex,'b','Linewidth',2)
hold on
plot(T,c2_qs,'r','Linewidth',2)
title('c2_{ex} (blue), c2_{u} (red)','Fontsize',24)


% function: eq_fun, exact dynamics
function dy=eq_fun(t,y,param)

e0=param(1);
k2=param(2);
k1=param(3);
km1=param(4);
k3=param(5);
km3=param(6);

s=y(1);
i=y(2);
c1=y(3);
c2=y(4);

e=e0-c1-c2;

dy=zeros(4,1);

dy(1)=-k1*s*e+km1*c1;
dy(2)=-k3*i*e+km3*c2;
dy(3)=k1*s*e-(km1+k2)*c1;
dy(4)=k3*i*e-km3*c2;

end


% function: eq_qs_fun, uniform approximation
function dy=eq_qs_fun(t,y,param)

e0=param(1);
k2=param(2);
k1=param(3);
km1=param(4);
k3=param(5);
km3=param(6);

s=y(1);
i=y(2);

K_i=km3/k3;
K_m=(km1+k2)/k1;

c1=e0*s./(s+K_m*(1+i/K_i));
c2=K_m/K_i*e0*i./(s+K_m*(1+i/K_i));

e=e0-c1-c2;

dy=zeros(2,1);

dy(1)=-k1*s*e+km1*c1;
dy(2)=-k3*i*e+km3*c2;

end
