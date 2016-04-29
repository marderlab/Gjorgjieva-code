cd ..
V1=dlmread('V1_mean_0_sig_1_gNa_227.052_gCaT_2.7781_gCaS_3.3469_gA_30.4321_gKCa_121.124_gK_75.649_gH_1.3194_gL_0.1631_gsyn_0.005_0.005_corr_0_tc_1_dt_0.01_num_2.dat');
V2=dlmread('V2_mean_0_sig_1_gNa_227.052_gCaT_2.7781_gCaS_3.3469_gA_30.4321_gKCa_121.124_gK_75.649_gH_1.3194_gL_0.1631_gsyn_0.005_0.005_corr_0_tc_1_dt_0.01_num_2.dat');

sp1=dlmread('spikes1_mean_0_sig_1_gNa_227.052_gCaT_2.7781_gCaS_3.3469_gA_30.4321_gKCa_121.124_gK_75.649_gH_1.3194_gL_0.1631_gsyn_0.005_0.005_corr_0_tc_1_dt_0.01_num_0.dat');
sp2=dlmread('spikes2_mean_0_sig_1_gNa_227.052_gCaT_2.7781_gCaS_3.3469_gA_30.4321_gKCa_121.124_gK_75.649_gH_1.3194_gL_0.1631_gsyn_0.005_0.005_corr_0_tc_1_dt_0.01_num_0.dat');

cd data_base

%%
cd ..
V1=dlmread('V1_mean_5_sig_1_gNa_227.052_gCaT_2.7781_gCaS_3.3469_gA_45_gKCa_121.124_gK_75.649_gH_1.3194_gL_0.1631_gsyn_0.2_0.2_corr_0_tc_1_dt_0.01_num_2.dat');
V2=dlmread('V2_mean_5_sig_1_gNa_227.052_gCaT_2.7781_gCaS_3.3469_gA_45_gKCa_121.124_gK_75.649_gH_1.3194_gL_0.1631_gsyn_0.2_0.2_corr_0_tc_1_dt_0.01_num_2.dat');

sp1=dlmread('spikes1_mean_5_sig_1_gNa_227.052_gCaT_2.7781_gCaS_3.3469_gA_45_gKCa_121.124_gK_75.649_gH_1.3194_gL_0.1631_gsyn_0.2_0.2_corr_0_tc_1_dt_0.01_num_0.dat');
sp2=dlmread('spikes2_mean_5_sig_1_gNa_227.052_gCaT_2.7781_gCaS_3.3469_gA_45_gKCa_121.124_gK_75.649_gH_1.3194_gL_0.1631_gsyn_0.2_0.2_corr_0_tc_1_dt_0.01_num_0.dat');

cd data_base
%%
% count = 0;
% for ii=5000:1000:60000
%     V11 = V1(650000-ii*10:end);
%     V12 = V2(650000-ii*10:end);
%     count = count+1;
%     chi(count) = burst_exclusion(V11,V12,0:0.1:ii);
% end
% 
% figure;plot(5000:1000:60000,chi,'o-')
% 
% 
% %try an average of multiple shorter ones
% count = 0;
% for ii=5000:100:64000
%     V11 = V1(ii*10:(ii+1000)*10);
%     V12 = V2(ii*10:(ii+1000)*10);
%     count = count+1;
%     chi2(count) = burst_exclusion(V11,V12,0:0.1:5000);
% end
% 
% figure;plot(5000:100:64000,chi2,'o-')
clear chi
tic
count = 0;
%remove first 5000 sec, transient period
%normalize so that spike trains start at 0
spp1 = sp1(sp1>5000)-5000; 
spp2 = sp2(sp2>5000)-5000;

%starting at 1000, compute the burst exclusion metric for increasing
%intervals of sizes: 1000, 2000, ... 2 000 000
for ii=5000:5000:2000000
    if mod(count,50)<0.1
        count
    end
    sp11 = spp1(spp1<ii);
    sp12 = spp2(spp2<ii);
    count = count+1;
    chi(count) = burst_exclusion(sp11,sp12,ii);
end
toc
figure;plot(5000:5000:2000000,chi,'o-')
%%
clear chi2

%remove first 5000 sec, transient period
%normalize so that spike trains start at 0
spp1 = sp1(sp1>5000)-5000; 
spp2 = sp2(sp2>5000)-5000;

%compute the burst exclusion metric for multiple shorter intervals of
%length 5000 which move with a sliding window of 100
count = 0;
for ii=5000:100:2000000-5000
    %normalize each time so that the new spike trains starts at 0 and ends
    %at 5000 (total duration of each short spike train interval
    sp11 = spp1(spp1>ii & spp1<ii+5000)-ii; 
    sp12 = spp2(spp2>ii & spp2<ii+5000)-ii;
    count = count+1;
    chi2(count) = burst_exclusion(sp11,sp12,5000);
end

figure;plot(5000:100:2000000-5000,chi2,'o-')

%the following is not good for extracting spikes, because voltage is
%recorded only every 0.1 ms
%
% tv1 = 0.1:0.1:65000;
% Vth = 15;
% %get cross-correlation
% %first get spike times
% dV1 = diff(V1);
% x1 = find(dV1(1:end-1)>0 & dV1(2:end)<0);   %find where derivative switches sign     
% x2 = find(V1>Vth); 
% x3 = intersect(x1+1,x2);
% sp1 = tv1(x3);
% 
% dV2 = diff(V2);
% x1 = find(dV2(1:end-1)>0 & dV2(2:end)<0);   %find where derivative switches sign     
% x2 = find(V2>Vth); 
% x3 = intersect(x1+1,x2);
% sp2 = tv1(x3);
%%

%remove first 5000 sec, transient period
%normalize so that spike trains start at 0
spp1 = sp1(sp1>5000)-5000; 
spp2 = sp2(sp2>5000)-5000;

% bin_width = 1; %ms, for binning spike trains to compute correlations
% offset = 0.5*10^4;
% time_vec = [0:bin_width:2000*10^3];
% T = 500; %s
% corr_vec = [-T:bin_width:T];
% cell1_binned = histc(sp1-offset,time_vec); %
% cell2_binned = histc(sp2-offset,time_vec); %
% cc = xcov(cell1_binned,cell2_binned,T/bin_width,'coef'); %
figure;subplot(131);%plot(corr_vec,cc);xlim([-T T]);%title(sprintf('tau=%g, corr=%g',tau,corr));
%
bin_width2 = 5; %ms, for binning spike trains to compute correlations
offset = 0.5*10^4;
time_vec = [0:bin_width2:2000*10^3];
T = 500; %s
corr_vec = [-T:bin_width2:T];
cell1_binned = histc(sp1-offset,time_vec); %
cell2_binned = histc(sp2-offset,time_vec); %
cc = xcov(cell1_binned,cell2_binned,T/bin_width2,'coef'); %
%hold on;plot(corr_vec,cc/(bin_width2/bin_width),'c');xlim([-T T])
plot(corr_vec,cc/bin_width2);xlim([-T T]);%title(sprintf('tau=%g, corr=%g',tau,corr));

cc1 = xcov(cell1_binned,cell1_binned,T/bin_width2,'coef'); %
subplot(132);plot(corr_vec,cc1,'r');xlim([-T T])

cc2 = xcov(cell2_binned,cell2_binned,T/bin_width2,'coef'); %
subplot(133);plot(corr_vec,cc2,'k');xlim([-T T])


%%
%indeed: the periods vary, sometimes synchrony, other times not
%maybe the longer the trial, there is convergence to a common value?
%but how long is long enough?

t_vec=0.1:0.1:65000;
figure;plot(t_vec,V1);hold on;plot(t_vec,V2);
for ii=1:65
    xlim([0 1000]+(ii-1)*1000)
    pause(1)
end