% PARAMETERS
kp = 50;
km = 0.1;
a_0 = 10;
b_0 = 5;

% RIGHT HAND SIDE (RHS) FUNCTION  
f = @(t,y) kp.*(a_0-y).*(b_0-y) - km*y;

% INITIAL CONDITION
y_0 = 0;

% TIME SPAN
t_0 = 0;
t_max = 1;
tt = [t_0 t_max];

% SOLVE
[T,Y] = ode23s(f,tt,y_0);

% PLOT
plot(T,Y) 

