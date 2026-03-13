% time interval
t0=0;
tmax=1;

% discretization parameter
h=0.001;

% enter paraters
k1=102;
km1=50;
k2=1;
k3=26;
km3=50;

% initial conditions 
s0=5;
e0=0.001;

% solution
options=odeset('RelTol',5.e-13 ,'AbsTol',1.e-13*ones(1,3),'InitialStep',1.e-5);
[T,Y]=ode15s(@eq_fun,[t0:h:tmax],[s0 0 0],options,e0,k1,km1,k2,k3,km3);

options=odeset('RelTol',5.e-13 ,'AbsTol',1.e-13,'InitialStep',1.e-5);
[T,Y_qs]=ode15s(@eq_qs_fun,[t0:h:tmax],s0,options,e0,k1,km1,k2,k3,km3);

s_ex=Y(:,1);
c1_ex=Y(:,2);
c2_ex=Y(:,3);

s_qs=Y_qs;

Km=(km1+k2)/k1;
Ki=km3/k3;

c1_qs=Ki*e0*s_qs./(Ki*Km+Ki*s_qs+s_qs.^2);
c2_qs=e0*s_qs.^2./(Ki*Km+Ki*s_qs+s_qs.^2);

% compute relative error between Y1 and Y2
err_s=norm(s_qs-s_ex,inf)/norm(s_ex,inf);
err_c1=norm(c1_qs-c1_ex,inf)/norm(c1_ex,inf);
err_c2=norm(c2_qs-c2_ex,inf)/norm(c2_ex,inf);
disp(['error on component s = ',num2str(err_s)]);
disp(['error on component c1 = ',num2str(err_c1)]);
disp(['error on component c2 = ',num2str(err_c2)]);

% plot
subplot(3,1,1)
plot(T,s_ex,'b','Linewidth',2)
hold on
plot(T,s_qs,'r','Linewidth',2)
title('s_{ex} (blue), s_{qs} (red)','Fontsize',24)

subplot(3,1,2)
plot(T,c1_ex,'b','Linewidth',2)
hold on
plot(T,c1_qs,'r','Linewidth',2)
title('c1_{ex} (blue), c1_{qs} (red)','Fontsize',24)

subplot(3,1,3)
plot(T,c2_ex,'b','Linewidth',2)
hold on
plot(T,c2_qs,'r','Linewidth',2)
title('c2_{ex} (blue), c2_{qs} (red)','Fontsize',24)


% function: eq_fun, exact dynamics
function dy=eq_fun(t,y,e0,k1,km1,k2,k3,km3)

s=y(1);
c1=y(2);
c2=y(3);

dy=zeros(3,1);

e=e0-c1-c2;

dy(1)=-k1*s*e+km1*c1-k3*s*c1+km3*c2;
dy(2)=k1*s*e-km1*c1-k2*c1-k3*s*c1+km3*c2;
dy(3)=k3*s*c1-km3*c2;

end


% function: eq_qs_fun, quasi-steady state approximation
function dy=eq_qs_fun(t,y,e0,k1,km1,k2,k3,km3)

s=y(1);

Km=(km1+k2)/k1;
Ki=km3/k3;

c1=Ki*e0*s./(Ki*Km+Ki*s+s.^2);
c2=e0*s.^2./(Ki*Km+Ki*s+s.^2);

e=e0-c1-c2;

dy=-k2*c1;

end
