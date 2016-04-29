cell1 = 894;%312;
cell2 = 1994;%1125;
g_syn1 = 0.06;%0.07;
g_syn2 = 0.03;%0.09;

Vsyn = -78;
tt = 3;
gmax = 100;
phi = 0.5;

% the tricky one... is this phase advance or delay??
% Vsyn = -78;
% tt = 1;
% gmax = 1;
% phi = 0.15;

% if (gmax > 0.99)
%     V_nosyn = dlmread(sprintf('PRC_data/V_in1_%d_in2_%d_gsyn_%g_%g_Vsyn_%g_gapp_%g_fracper_%g_phase_%g_cell_%d.dat',cell1,cell2,g_syn1,g_syn2,Vsyn,gmax,Tdur_vec(tt),phi,0));
% else
%     V_nosyn = dlmread(sprintf('PRC_data/V_in1_%d_in2_%d_gsyn_%g_%g_Vsyn_%g_gapp_%g_fracper_%g_phase_%g_cell_%d.dat',cell1,cell2,g_syn1,g_syn2,Vsyn,gmax,Tdur_vec(tt),phi,2));
% end

%input into cell 1
%V_nosyn = dlmread(sprintf('PRC_data/V_in1_%d_in2_%d_gsyn_%g_%g_Vsyn_%g_gapp_%g_fracper_%g_phase_%g_cell_%d.dat',cell1,cell2,g_syn1,g_syn2,Vsyn,gmax,Tdur_vec(tt),phi,0));

%input into cell 2
%V_nosyn = dlmread(sprintf('PRC_data/V_in1_%d_in2_%d_gsyn_%g_%g_Vsyn_%g_gapp_%g_fracper_%g_phase_%g_cell_%d.dat',cell1,cell2,g_syn1,g_syn2,Vsyn,gmax,Tdur_vec(tt),phi,1));

%input into cells 1&2
V_nosyn = dlmread(sprintf('PRC_data/V_in1_%d_in2_%d_gsyn_%g_%g_Vsyn_%g_gapp_%g_fracper_%g_phase_%g_cell_%d.dat',cell1,cell2,g_syn1,g_syn2,Vsyn,gmax,Tdur_vec(tt),phi,2));



V1 = V_nosyn(V_nosyn(:,1)==0,2);
dV1 = diff(V1); 
x1 = find(dV1(1:end-1)>0 & dV1(2:end)<0);   %find where derivative switches sign     
x2 = find(V1>Vth); 
x3 = intersect(x1+1,x2);
sp1 = tot_T(x3); %get spike times
[CV1, num_sp1, period1, burst_dur1, first_spike1, last_spike1] = get_train_properties_HCO(sp1,0,30);
time_burst_start_ind1 = find(first_spike1>2000,2);
(period1(time_burst_start_ind1(1)) - ref_period1)/ref_period1

figure;plot(0:0.1:5000,V1,'linewidth',2);hold on;plot(0:0.1:5000,Vs1,'--','linewidth',2);xlim([1900 2600]) %Vs1 or V0_1
set(gcf,'position',[0 0 1200 200]);
set(gca,'fontsize',20)
ylim([-80 60])

V2 = V_nosyn(V_nosyn(:,1)==1,2);
dV2 = diff(V2); 
x1 = find(dV2(1:end-1)>0 & dV2(2:end)<0);   %find where derivative switches sign     
x2 = find(V2>Vth); 
x3 = intersect(x1+1,x2);
sp2 = tot_T(x3); %get spike times
[CV2, num_sp2, period2, burst_dur2, first_spike2, last_spike2] = get_train_properties_HCO(sp2,0,30);
time_burst_start_ind2 = find(first_spike2>2000,2);
(period2(time_burst_start_ind2(1)) - ref_period2)/ref_period2

figure;plot(0:0.1:5000,V2,'linewidth',2);hold on;plot(0:0.1:5000,Vs2,'--','linewidth',2);xlim([1900 2600]) %Vs2 or V0_2
set(gcf,'position',[0 0 1200 200]);
set(gca,'fontsize',20)
ylim([-80 60])


%saveas(gcf,sprintf('V_in1_%d_in2_%d_gsyn_%g_%g_Vsyn_%g_gapp_%g_fracper_%g_phase_%g_cell_%d.fig',cell1,cell2,g_syn1,g_syn2,Vsyn,gmax,Tdur_vec(tt),phi,0),'fig');