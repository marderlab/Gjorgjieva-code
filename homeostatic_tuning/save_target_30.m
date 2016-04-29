
dt = 0.1;
taug = 1.0e2;
taum = 10*taug;
tstop = 1200*taug;
ms_per_samp = tstop/2000;
t=dt:ms_per_samp:tstop;

load plotcolours;
colrs = brighten(colrs,0.7);
linw = 1.0;
figure;
g1_ave = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mean = 5;
sig = 5;
tau = 500;
coef = 0.5;
Ntot = 50;
numm = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%for ii=setdiff([1:Ntot],[21]) %for tau=100 ms
%for ii=setdiff([1:Ntot],[6 14 32 36]) %for tau=200 ms
for ii=setdiff([1:Ntot],[2 6 10 30 32 33 34 35 37 39 44 46]) %for tau=500 ms
    %set(gcf,'Renderer','Zbuffer');
    %ii
    hold on;
    set(gca,'colororder',colrs);
    g1=dlmread(sprintf('g1_mean_%g_sig_%g_tau_%g_coef_%g_target_30_taug_100_taum_1000_num_%d.dat',mean,sig,tau,coef,ii));
    %g1=dlmread(sprintf('g1_target_30_taug_100_taum_1000_num_%d.dat',ii));
    %g2=dlmread(sprintf('g2_target_30_taug_100_taum_1000_num_%d.dat',ii));
    t = linspace(1/taug,tstop/taug,size(g1,1));
    g1(1,9) = g1(1,8); %for the IC, I omitted to save g_leak but I saved g_syn in the place of g_leak
    g1 = g1(:,[1:7 9]);
    plot(t,g1,'linewidth',linw);
    %set(gca,'yscale','log','xscale','log');
    box off;
    axis tight;
    drawnow;
    g1_ave = g1_ave+g1;
    numm = numm + 1;
    %pause(1)
end
load plotcolours;
linw = 2.0;
%set(gcf,'Renderer','Zbuffer');
set(gca,'colororder',colrs);
plot(t,g1_ave/numm,'linewidth',linw);
%set(gca,'ytick',logspace(-3,3,3),'xminortick','on','yminortick','on');
xlabel('time/tau_g');
ylabel('conductance (\muS)');
tidyfonts(19);
sizefig(700,400);


