
V_pre = [-100:0.1:50];

V_half = -45;
V_slope = 5;
tau_syn = 100;
S_inf2 = tanh((V_pre-V_half)./V_slope);
S_inf = S_inf2.*(V_pre>V_half);

figure;plot(V_pre,S_inf,'b','linewidth',2);hold on;%plot(V_pre,S_inf2,'r--')

V_half = 40;
V_slope = -2;
S_inf_new = 1.0./(1.0+exp((V_pre + V_half)/V_slope)); 
plot(V_pre,S_inf_new,'k--','linewidth',2)

%%
G = 10;
V_post = [-100:0.1:50];
E_syn = -70;
S = rand(1,length(V_post));%linspace(0,1,length(V_post));%
I_syn = G.*S.*(E_syn-V_post);

figure;plot(V_post,I_syn)
%%
clear Svar
N = 10^3;

V1 = dlmread('dat_files/V1_mean_0_sig_0_gNa_227.052_gCaT_2.7781_gCaS_3.3469_gA_30.4321_gKCa_121.124_gK_75.649_gH_1.3194_gL_0.1631_gsyn_0.001_0.001_corr_0_tc_1_dt_0.01_num_2.dat');
V2 = dlmread('dat_files/V2_mean_0_sig_0_gNa_227.052_gCaT_2.7781_gCaS_3.3469_gA_30.4321_gKCa_121.124_gK_75.649_gH_1.3194_gL_0.1631_gsyn_0.001_0.001_corr_0_tc_1_dt_0.01_num_2.dat');

Svar(1) = 0;
Svar2(1) = 0;
dt = 1;
for nn=2:N
     S_inf = tanh((V1(nn)-V_half)/V_slope).*(V1(nn)>V_half);
     Svar(nn) = Svar(nn-1) + (1-exp(-dt/(tau_syn*(1-S_inf))))*(S_inf - Svar(nn-1));
     Svar2(nn) = Svar2(nn-1) + (1-exp(-dt/(tau_syn)))*(S_inf - Svar2(nn-1));
end
I_syn = G.*Svar.*(E_syn-V2(1:N)');


figure;subplot(3,1,1);plot(V1(1:N),'b');
subplot(3,1,2);plot(Svar,'r');hold on;plot(Svar2,'k')
subplot(3,1,3);plot(I_syn,'k');