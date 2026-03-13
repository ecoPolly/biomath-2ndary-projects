% time interval
t0=0;
tmax=100;

% discretization parameter
h=0.001;

% enter paraters
k1=4e+6;
km1=25;
k2=15;
K_m=1e-5;

% initial conditions 
s0=K_m;
e0_v=[1e-1 1e-2 1e-3 1e-4]*K_m;

% transient time constant
tau=1/(k1*s0);

% solution
for i=1:4
   e0=e0_v(i);
   disp(['e0 = ',num2str(e0)]);

   options=odeset('RelTol',5.e-13 ,'AbsTol',1.e-13*ones(1,2),'InitialStep',1.e-5);
   [T,Y]=ode15s(@eq_fun,[t0:h:tmax],[s0 0],options,e0,k1,km1,k2,K_m);

   options=odeset('RelTol',5.e-13 ,'AbsTol',1.e-13,'InitialStep',1.e-5);
   [T,Y_u]=ode15s(@eq_u_fun,[t0:h:tmax],s0,options,e0,k1,km1,k2,K_m);

   s_ex=Y(:,1);
   c_ex=Y(:,2);

   s_u=Y_u;
   c_u=e0*s_u./(K_m+s_u)-e0*s0./(K_m+s0).*exp(-(K_m+s0)*k1*T);

   % compute relative error between Y and Y_u
   err_s=norm(s_u-s_ex,inf)/norm(s_ex,inf);
   err_c=norm(c_u-c_ex,inf)/norm(c_ex,inf);
   disp(['error on component s = ',num2str(err_s)]);
   disp(['error on component c = ',num2str(err_c)]);

   % plot
   figure
   subplot(2,2,1)
   plot(T,s_ex,'b','Linewidth',2)
   hold on
   plot(T,s_u,'r','Linewidth',2)
   title('s_{ex} (blue), s_{u} (red)','Fontsize',24)

   subplot(2,2,2)
   plot(T,s_ex,'b','Linewidth',2)
   hold on
   plot(T,s_u,'r','Linewidth',2)
   axis([0 2*tau min(s_ex) max(s_ex)])
   title('s_{ex} (blue), s_{u} (red)','Fontsize',24)


   subplot(2,2,3)
   plot(T,c_ex,'b','Linewidth',2)
   hold on
   plot(T,c_u,'r','Linewidth',2)
   title('c_{ex} (blue), c_{u} (red)','Fontsize',24)

   subplot(2,2,4)
   plot(T,c_ex,'b','Linewidth',2)
   hold on
   plot(T,c_u,'r','Linewidth',2)
   axis([0 2*tau min(c_ex) max(c_ex)])
   title('c_{ex} (blue), c_{u} (red)','Fontsize',24)
end

% function: eq_fun, exact dynamics
function dy=eq_fun(t,y,e0,k1,km1,k2,K_m)

s=y(1);
c=y(2);

dy=zeros(2,1);

e=e0-c;

dy(1)=-k1*s*e+km1*c;
dy(2)=k1*s*e-(km1+k2)*c;

end


% function: eq_u_fun, uniform approximation
function dy=eq_u_fun(t,y,e0,k1,km1,k2,K_m)

s=y(1);

dy=-k2*e0*s./(K_m+s);

end
