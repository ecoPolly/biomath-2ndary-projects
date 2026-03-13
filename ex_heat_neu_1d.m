% INPUT 
n = input('enter the number of space intervals = ');
m = input('enter the number of time intervals = ');

% PARAMETERS
L=1;
T=0.8;
h=L/n;
k=T/m;
x=linspace(0,L,n+1)';
tt=linspace(0,T,m+1)';

% ASSEMBLE MATRIX 
AA=2*diag(ones(n+1,1))-diag(ones(n,1),1)-diag(ones(n,1),-1);
AA(1,1)=1;
AA(end,end)=1;
AA=1/h^2*AA;
MM=eye(n+1);
MM(1,1)=1/2;
MM(end,end)=1/2;
AA=1/k*MM+AA;

% ASSEMBLE IC
u0_fun=@(x) 100*cos(pi*x);  
uu0=u0_fun(x);

% SET EXACT SOLUTION
sol_ex=@(x,t) 100*exp(-pi*pi*t)*cos(pi*x);

% SOLVE
uu=zeros(n+1,1);
uex=sol_ex(x,tt(1));
err(1)=norm(uu0-uex,inf)/norm(uex,inf);

for j=1:m-1
   bb=1/k*MM*uu0;
   uu=AA\bb;

   plot(x,uu,'r-o')
   hold off
   axis([0 L -110 110])
   pause(0.1);

   uex=sol_ex(x,tt(j+1));
   err(j+1)=norm(uu-uex,inf)/norm(uex,inf);

   uu0=uu;
end
 
% PRINT ERROR
err=max(err)

 
