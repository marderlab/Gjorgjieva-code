sp1 = dlmread('spikes1_mean_0_sig_0_gNa_227.052_gCaT_2.7781_gCaS_3.3469_gA_30.4321_gKCa_121.124_gK_75.649_gH_1.3194_gL_0.1631_gsyn_0.005_0.005_corr_0_tc_1_dt_0.01_num_2.dat');
sp2 = dlmread('spikes2_mean_0_sig_0_gNa_227.052_gCaT_2.7781_gCaS_3.3469_gA_30.4321_gKCa_121.124_gK_75.649_gH_1.3194_gL_0.1631_gsyn_0.005_0.005_corr_0_tc_1_dt_0.01_num_2.dat');
%sp1 = dlmread('spikes1_mean_0_sig_0_gNa_227.052_gCaT_2.2781_gCaS_2.8469_gA_30.4321_gKCa_1121.12_gK_75.649_gH_1.3194_gL_0.1631_gsyn_0.01_0.01_corr_0_tc_1_dt_0.01_num_2.dat');
%sp2 = dlmread('spikes2_mean_0_sig_0_gNa_227.052_gCaT_2.2781_gCaS_2.8469_gA_30.4321_gKCa_1121.12_gK_75.649_gH_1.3194_gL_0.1631_gsyn_0.01_0.01_corr_0_tc_1_dt_0.01_num_2.dat');

bin_width = 0.01; %s, for binning spike trains to compute correlations
offset = 0.5*10^4;

time_vec = [0:bin_width:6*10^4];
T = 150; %s
corr_vec = [-T:bin_width:T];
%sp1 = sp1(sp1>offset); sp2 = sp2(sp2>offset);

cell1_binned = histc(sp1-offset,time_vec); %
cell2_binned = histc(sp2-offset,time_vec); %
%cell2_binned = histc((sp1+10)-offset+0*rand(length(sp1),1),time_vec); %
cc = xcov(cell1_binned,cell2_binned,T/bin_width); %
figure;plot(corr_vec,cc);xlim([-T T])


cc1 = xcov(cell1_binned,cell1_binned,T/bin_width); %
figure;plot(corr_vec,cc1);xlim([-T T])

cc2 = xcov(cell2_binned,cell2_binned,T/bin_width); %
figure;plot(corr_vec,cc2);xlim([-T T])

%%
clear spp
jj = 0;
dt = 0.01;
rate = 2;
for ii=1:6*10^4
    r = rand;
    if (rate*dt > rand)
        jj = jj+1;
        spp(jj) = ii*dt;
    end
end

cell1_binned = histc(spp,time_vec); %
cell2_binned = histc((spp+10),time_vec); %
cc1 = xcov(cell1_binned,cell1_binned,2000); %
figure;plot(cc1)
cc2 = xcov(cell1_binned,cell2_binned,2000); %
figure;plot(cc2)