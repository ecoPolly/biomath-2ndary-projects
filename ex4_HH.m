% EX 4

% Initial conditions (resting)
y0=[2.7570e-04   5.2934e-02   5.9611e-01   3.1768e-01]';

T=100;
Tsti=200;
di=200;

Iapp=[2.5 3 4 5 6 7 8 9 10 15 20 25 50 100 200];

for i=1:length(Iapp)
    [t,y]=ode15s(@HH,[0 T],y0,[],di,Tsti,Iapp(i));
    figure(i)
    plot(t,y(:,1))
    hold on
end

