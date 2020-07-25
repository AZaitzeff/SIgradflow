addpath('../../../MATLAB/PlotPub/lib/')
[trueU,~,~] = initializebarenblatt(N,2+1/4);
%[U,x,h]=initializebump(N);
[U,x,h] = initializebarenblatt(N,1/4);

plot(x,U,x,trueU);hold on
load('pmewm2o16385N32nt.mat');
plot(x,U)
load('pmewm3o16385N32nt.mat');
plot(x,U)

plt=Plot();
plt.Colors = {[.5,.5,.5],[0,0,0],[.5,.5,.5],[0,0,0]};
plt.LineStyle = {'-','-','--',':'};
plt.XLabel='x';
plt.YLabel='y';
plt.Legend = {'Initial', 'Final', '2nd Order','3rd Order'};
plt.XLim=[-3,3];
plt.export('barenblatt.png');
close


load('pmewm2o16385N32nt.mat');
ts=0:dt:T;
plot(ts,energy);hold on
load('pmewm3o16385N32nt.mat');
plot(ts,energy)
plt=Plot();
plt.Colors = {[.5,.5,.5],[0,0,0]};
plt.LineStyle = {'-','--'};
plt.XLabel='t';
plt.YLabel='E(u)';
plt.Legend = {'2nd Order','3rd Order'};
plt.export('barenblattenergy.png');
close