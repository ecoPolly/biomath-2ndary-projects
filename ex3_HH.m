% EX 3

% Initial conditions (resting)
y0=[2.7570e-04   5.2934e-02   5.9611e-01   3.1768e-01]';

T=100;
Iapp=100;
Tsti=0.1;

di=[20 17 15 12 10 8 7 6 5 4 3 2 1.5];

for i=1:length(di)
  tt=0; 
  yy=y0;
  while(tt<=T)
    [t,y]=ode15s(@HH,[tt tt+di(i)],yy,[],di(i),Tsti,Iapp);
    figure(i)
    plot(t,y(:,1))
    hold on
    tt=t(end);
    yy=y(end,:);
  end
end
