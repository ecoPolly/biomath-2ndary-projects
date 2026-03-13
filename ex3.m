% PARAMETERS
r_p = 0.6;
r_s = -0.2;
sigma_p = 2e-6;
sigma_s = 1e-4;
alpha = 1e-3;
mu = 0.01;

% INITIAL CONDITION
y_0 = [5000 40];

% TIME SPAN
t_0 = 0;
t_max = 80;
tt = [t_0 t_max];

% CASE 1
b_d = 0;
param = [r_p r_s sigma_p sigma_s alpha mu b_d];

% SOLVE
[T,Y] = ode45(@fun_sys,tt,y_0,[],param);

% PLOT
subplot(1,2,1)
plot(T,Y(:,1),'b','Linewidth',3) 
hold on
subplot(1,2,2)
plot(T,Y(:,2),'b','Linewidth',3)
hold on

% CASE 2
b_d = 0.1;
param = [r_p r_s sigma_p sigma_s alpha mu b_d];

% SOLVE
[T,Y] = ode45(@fun_sys,tt,y_0,[],param);

% PLOT
subplot(1,2,1)
plot(T,Y(:,1),'r','Linewidth',3)
hold on
subplot(1,2,2)
plot(T,Y(:,2),'r','Linewidth',3)
hold on

% CASE 3
b_d = 0.3;
param = [r_p r_s sigma_p sigma_s alpha mu b_d];

% SOLVE
[T,Y] = ode45(@fun_sys,tt,y_0,[],param);

% PLOT
subplot(1,2,1)
plot(T,Y(:,1),'g','Linewidth',3)
hold on
subplot(1,2,2)
plot(T,Y(:,2),'g','Linewidth',3)
hold on


function dy = fun_sys(t,y,param)

% PARAMETERS
r_p = param(1);
r_s = param(2);
sigma_p = param(3);
sigma_s = param(4);
alpha = param(5);
mu = param(6);
b_d = param(7);

% EXTRACT VARIABLES
p = y(1);
s = y(2);

% SET RHS
dy = zeros(2,1);

dy(1) = (r_p -sigma_p*p - alpha*s -b_d)*p;
dy(2) = (r_s -sigma_s*s + mu*alpha*p -b_d)*s;

end  
