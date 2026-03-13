% PARAMETERS
r = 0.06;
K = 150;
b = @(t) -0.03*(2-exp(-5*sin(pi*t).^2));

% RIGHT HAND SIDE (RHS) FUNCTION  
f = @(t,y) r*y-r/K*y.^2+b(t).*y;

% INITIAL CONDITION
y_0 = 100;

% TIME SPAN
t_0 = 2017.5;
t_max = 2025.5;
h = 1/12;
tt = [t_0:h:t_max];

% SET OPTIONS
options = odeset('RelTol',5.e-13 ,'AbsTol',[1.e-13],'InitialStep',1.e-5,'MaxStep',5);

% SOLVE
[T,Y] = ode45(f,tt,y_0);

% PLOT
plot(T,Y) 
hold on

% EE
N = length(tt);
uu = zeros(N,1);
uu(1) = y_0;

for i = 1:N-1
   uu(i+1) = uu(i)+h*f(tt(i),uu(i));
end
plot(tt,uu,'r')
uu(end)

% EE
h = 0.1/12;
tt = [t_0:h:t_max];

N = length(tt);
uu = zeros(N,1);
uu(1) = y_0;

for i = 1:N-1
   uu(i+1) = uu(i)+h*f(tt(i),uu(i));
end
plot(tt,uu,'g')
uu(end)


% FUNCTIONS
function dy = fun_sys(t,y)

% PARAMETERS
r = 0.06;
K = 150;
bt = -0.03*(2-exp(-5*sin(pi*t).^2));

dy = r*y-r/K*y.^2+bt.*y;

end
