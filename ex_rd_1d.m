% INPUT 
n = input('enter the number of intervals = ');

% PARAMETERS
f=@(x) (pi^2+1)*cos(pi*x);
h=1/n;
x=linspace(0,1,n+1)';

% ASSEMBLE MATRIX 
AA=2*diag(ones(n+1,1))-diag(ones(n,1),1)-diag(ones(n,1),-1);
AA(1,1)=1;
AA(end,end)=1;
AA=1/h^2*AA;
MM=eye(n+1);
MM(1,1)=1/2;
MM(end,end)=1/2;
AA=AA+MM;

% ASSEMBLE RHS
bb=f(x);
bb=MM*bb;

% SOLVE
uu=zeros(n+1,1);
uu=AA\bb;

% ERROR
sol_ex=@(x) cos(pi*x);
uex=sol_ex(x);
err=norm(uu-uex,inf)/norm(uex,inf)

% PLOT
plot(x,uex,x,uu,'r')
 
