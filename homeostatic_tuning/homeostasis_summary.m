%So when DC=0, need to use a small-ish target and small tau_syni which
%results in small synaptic conductances. (Beacuse ratio of time constants 
%<=>ratio of conductances). Indeed, activity is bursty but not
%always alternating suggesting the lack of robustness of this model. It is
%intrisically dominated because removing synaptic connections in the tuned
%models still preserves bursting activity.

%But when DC=5 and using large-ish target and large tau_syni which results
%in larger synaptic conductances, can generate very nice bursts but
%intrinsically the neurons are spike-y. 

%Still need to do: tuning in the presense of noisy input


%DC=0, target 10, tau_syn 0.001
%1.
V1=dlmread('V1_mean_0_sig_0_gNa_206.56_gCaT_2.63618_gCaS_3.17481_gA_27.6294_gKCa_110.184_gK_68.9003_gH_1.14958_gL_0.072641_gsyn_0.014247_0.017793_corr_0_tc_0_dt_0.01_num_1.dat');
V2=dlmread('V2_mean_0_sig_0_gNa_257.751_gCaT_3.28653_gCaS_3.92129_gA_34.5_gKCa_137.44_gK_85.9744_gH_1.49074_gL_0.158543_gsyn_0.014247_0.017793_corr_0_tc_0_dt_0.01_num_1.dat');
figure;plot(0:0.1:5000,V1);hold on;plot(0:0.1:5000,V2)
xlim([4000 4800])
%without synapse
V1=dlmread('V1_mean_0_sig_0_gNa_206.56_gCaT_2.63618_gCaS_3.17481_gA_27.6294_gKCa_110.184_gK_68.9003_gH_1.14958_gL_0.072641_gsyn_0_0_corr_0_tc_0_dt_0.01_num_1.dat');
V2=dlmread('V2_mean_0_sig_0_gNa_257.751_gCaT_3.28653_gCaS_3.92129_gA_34.5_gKCa_137.44_gK_85.9744_gH_1.49074_gL_0.158543_gsyn_0_0_corr_0_tc_0_dt_0.01_num_1.dat');
figure;plot(0:0.1:5000,V1);hold on;plot(0:0.1:5000,V2)
xlim([4000 4800])

%2.
V1=dlmread('V1_mean_0_sig_0_gNa_261.599_gCaT_3.17728_gCaS_3.94158_gA_34.9997_gKCa_139.497_gK_87.2069_gH_1.54282_gL_0.1492_gsyn_0.017946_0.018085_corr_0_tc_0_dt_0.01_num_1.dat');
V2=dlmread('V2_mean_0_sig_0_gNa_259.265_gCaT_3.24277_gCaS_3.80731_gA_34.5619_gKCa_138.253_gK_86.5174_gH_1.44996_gL_0.154454_gsyn_0.017946_0.018085_corr_0_tc_0_dt_0.01_num_1.dat');
figure;plot(0:0.1:5000,V1);hold on;plot(0:0.1:5000,V2)
xlim([4000 4800])
%without synapse
V1=dlmread('V1_mean_0_sig_0_gNa_261.599_gCaT_3.17728_gCaS_3.94158_gA_34.9997_gKCa_139.497_gK_87.2069_gH_1.54282_gL_0.1492_gsyn_0_0_corr_0_tc_0_dt_0.01_num_1.dat');
V2=dlmread('V2_mean_0_sig_0_gNa_259.265_gCaT_3.24277_gCaS_3.80731_gA_34.5619_gKCa_138.253_gK_86.5174_gH_1.44996_gL_0.154454_gsyn_0_0_corr_0_tc_0_dt_0.01_num_1.dat');
figure;plot(0:0.1:5000,V1);hold on;plot(0:0.1:5000,V2)
xlim([4000 4800])

%3.
V1=dlmread('V1_mean_0_sig_0_gNa_224.324_gCaT_2.84044_gCaS_3.30154_gA_29.9619_gKCa_119.686_gK_74.7952_gH_1.38842_gL_0.09902_gsyn_0.015832_0.018278_corr_0_tc_0_dt_0.01_num_1.dat');
V2=dlmread('V2_mean_0_sig_0_gNa_270.157_gCaT_3.29957_gCaS_3.9953_gA_36.0677_gKCa_144.064_gK_90.1214_gH_1.58027_gL_0.172585_gsyn_0.015832_0.018278_corr_0_tc_0_dt_0.01_num_1.dat');
figure;plot(0:0.1:5000,V1);hold on;plot(0:0.1:5000,V2)
xlim([4000 4800])
%without synapse
V1=dlmread('V1_mean_0_sig_0_gNa_224.324_gCaT_2.84044_gCaS_3.30154_gA_29.9619_gKCa_119.686_gK_74.7952_gH_1.38842_gL_0.09902_gsyn_0_0_corr_0_tc_0_dt_0.01_num_1.dat');
V2=dlmread('V2_mean_0_sig_0_gNa_270.157_gCaT_3.29957_gCaS_3.9953_gA_36.0677_gKCa_144.064_gK_90.1214_gH_1.58027_gL_0.172585_gsyn_0_0_corr_0_tc_0_dt_0.01_num_1.dat');
figure;plot(0:0.1:5000,V1);hold on;plot(0:0.1:5000,V2)
xlim([4000 4800])

%4.
V1=dlmread('V1_mean_0_sig_0_gNa_267.968_gCaT_3.21155_gCaS_4.10866_gA_35.7187_gKCa_142.901_gK_89.3904_gH_1.53632_gL_0.180901_gsyn_0.018743_0.018088_corr_0_tc_0_dt_0.01_num_1.dat');
V2=dlmread('V2_mean_0_sig_0_gNa_260.162_gCaT_3.25607_gCaS_3.9921_gA_34.8482_gKCa_138.71_gK_86.8409_gH_1.51444_gL_0.13436_gsyn_0.018743_0.018088_corr_0_tc_0_dt_0.01_num_1.dat');
figure;plot(0:0.1:5000,V1);hold on;plot(0:0.1:5000,V2)
xlim([4000 4800])
%without synapse
V1=dlmread('V1_mean_0_sig_0_gNa_267.968_gCaT_3.21155_gCaS_4.10866_gA_35.7187_gKCa_142.901_gK_89.3904_gH_1.53632_gL_0.180901_gsyn_0_0_corr_0_tc_0_dt_0.01_num_1.dat');
V2=dlmread('V2_mean_0_sig_0_gNa_260.162_gCaT_3.25607_gCaS_3.9921_gA_34.8482_gKCa_138.71_gK_86.8409_gH_1.51444_gL_0.13436_gsyn_0_0_corr_0_tc_0_dt_0.01_num_1.dat');
figure;plot(0:0.1:5000,V1);hold on;plot(0:0.1:5000,V2)
xlim([4000 4800])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DC=0, target 12, tau_syn 0.001
%1.
V1=dlmread('V1_mean_0_sig_0_gNa_209.19_gCaT_2.64534_gCaS_3.0919_gA_28.0673_gKCa_111.58_gK_69.788_gH_1.15657_gL_0.001979_gsyn_0.014188_0.018444_corr_0_tc_0_dt_0.01_num_1.dat');
V2=dlmread('V2_mean_0_sig_0_gNa_266.806_gCaT_3.22875_gCaS_3.93999_gA_35.5702_gKCa_142.371_gK_88.9916_gH_1.48178_gL_0.110801_gsyn_0.014188_0.018444_corr_0_tc_0_dt_0.01_num_1.dat');
figure;plot(0:0.1:5000,V1);hold on;plot(0:0.1:5000,V2)
xlim([4000 4800])
%without synapse
V1=dlmread('V1_mean_0_sig_0_gNa_209.19_gCaT_2.64534_gCaS_3.0919_gA_28.0673_gKCa_111.58_gK_69.788_gH_1.15657_gL_0.001979_gsyn_0_0_corr_0_tc_0_dt_0.01_num_1.dat');
V2=dlmread('V2_mean_0_sig_0_gNa_266.806_gCaT_3.22875_gCaS_3.93999_gA_35.5702_gKCa_142.371_gK_88.9916_gH_1.48178_gL_0.110801_gsyn_0_0_corr_0_tc_0_dt_0.01_num_1.dat');
figure;plot(0:0.1:5000,V1);hold on;plot(0:0.1:5000,V2)
xlim([4000 4800])

%2.
V1=dlmread('V1_mean_0_sig_0_gNa_210.78_gCaT_2.6659_gCaS_3.12656_gA_28.1448_gKCa_112.419_gK_70.3823_gH_1.28029_gL_0.011775_gsyn_0.014967_0.014673_corr_0_tc_0_dt_0.01_num_1.dat');
V2=dlmread('V2_mean_0_sig_0_gNa_211.966_gCaT_2.60111_gCaS_3.16942_gA_28.2757_gKCa_113.044_gK_70.6234_gH_1.21927_gL_0.003883_gsyn_0.014967_0.014673_corr_0_tc_0_dt_0.01_num_1.dat');
figure;plot(0:0.1:5000,V1);hold on;plot(0:0.1:5000,V2)
xlim([4000 4800])
%without synapse
V1=dlmread('V1_mean_0_sig_0_gNa_210.78_gCaT_2.6659_gCaS_3.12656_gA_28.1448_gKCa_112.419_gK_70.3823_gH_1.28029_gL_0.011775_gsyn_0_0_corr_0_tc_0_dt_0.01_num_1.dat');
V2=dlmread('V2_mean_0_sig_0_gNa_211.966_gCaT_2.60111_gCaS_3.16942_gA_28.2757_gKCa_113.044_gK_70.6234_gH_1.21927_gL_0.003883_gsyn_0_0_corr_0_tc_0_dt_0.01_num_1.dat');
figure;plot(0:0.1:5000,V1);hold on;plot(0:0.1:5000,V2)
xlim([4000 4800])

%3.
V1=dlmread('V1_mean_0_sig_0_gNa_231.166_gCaT_2.87938_gCaS_3.4159_gA_30.9514_gKCa_123.319_gK_77.1471_gH_1.38725_gL_0.035871_gsyn_0.016302_0.020323_corr_0_tc_0_dt_0.01_num_1.dat');
V2=dlmread('V2_mean_0_sig_0_gNa_291.681_gCaT_3.59208_gCaS_4.47219_gA_39.0468_gKCa_155.582_gK_97.2542_gH_1.6_gL_0.15788_gsyn_0.016302_0.020323_corr_0_tc_0_dt_0.01_num_1.dat');
figure;plot(0:0.1:5000,V1);hold on;plot(0:0.1:5000,V2)
xlim([4000 4800])
%without synapse
V1=dlmread('V1_mean_0_sig_0_gNa_231.166_gCaT_2.87938_gCaS_3.4159_gA_30.9514_gKCa_123.319_gK_77.1471_gH_1.38725_gL_0.035871_gsyn_0_0_corr_0_tc_0_dt_0.01_num_1.dat');
V2=dlmread('V2_mean_0_sig_0_gNa_291.681_gCaT_3.59208_gCaS_4.47219_gA_39.0468_gKCa_155.582_gK_97.2542_gH_1.6_gL_0.15788_gsyn_0_0_corr_0_tc_0_dt_0.01_num_1.dat');
figure;plot(0:0.1:5000,V1);hold on;plot(0:0.1:5000,V2)
xlim([4000 4800])

%4. 
V1=dlmread('V1_mean_0_sig_0_gNa_257.055_gCaT_3.11332_gCaS_3.94067_gA_34.3175_gKCa_137.014_gK_85.7478_gH_1.50907_gL_0.097864_gsyn_0.017273_0.021_corr_0_tc_0_dt_0.01_num_1.dat');
V2=dlmread('V2_mean_0_sig_0_gNa_312.046_gCaT_3.75016_gCaS_4.69264_gA_41.6851_gKCa_166.45_gK_104.063_gH_1.78661_gL_0.189671_gsyn_0.017273_0.021_corr_0_tc_0_dt_0.01_num_1.dat');
figure;plot(0:0.1:5000,V1);hold on;plot(0:0.1:5000,V2)
xlim([4000 4800])
%without synapse
title('DC=0; tau_isyn=0.001, target=12, RUN 4')
V1=dlmread('V1_mean_0_sig_0_gNa_257.055_gCaT_3.11332_gCaS_3.94067_gA_34.3175_gKCa_137.014_gK_85.7478_gH_1.50907_gL_0.097864_gsyn_0_0_corr_0_tc_0_dt_0.01_num_1.dat');
V2=dlmread('V2_mean_0_sig_0_gNa_312.046_gCaT_3.75016_gCaS_4.69264_gA_41.6851_gKCa_166.45_gK_104.063_gH_1.78661_gL_0.189671_gsyn_0_0_corr_0_tc_0_dt_0.01_num_1.dat');
figure;plot(0:0.1:5000,V1);hold on;plot(0:0.1:5000,V2)
xlim([4000 4800])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DC=5; target=30; tau_syn=0.01
%1.
V1=dlmread('V1_mean_5_sig_0_gNa_852.746_gCaT_10.333_gCaS_12.5517_gA_113.63_gKCa_454.734_gK_284.132_gH_4.70886_gL_0.117436_gsyn_0.568939_0.561596_corr_0_tc_0_dt_0.01_num_1.dat');
V2=dlmread('V2_mean_5_sig_0_gNa_842.98_gCaT_10.1081_gCaS_12.5277_gA_112.336_gKCa_449.574_gK_281.021_gH_4.66558_gL_0.061134_gsyn_0.568939_0.561596_corr_0_tc_0_dt_0.01_num_1.dat');
figure;plot(0:0.1:5000,V1);hold on;plot(0:0.1:5000,V2)
xlim([4000 4500])
%without synapse
V1=dlmread('V1_mean_5_sig_0_gNa_852.746_gCaT_10.333_gCaS_12.5517_gA_113.63_gKCa_454.734_gK_284.132_gH_4.70886_gL_0.117436_gsyn_0_0_corr_0_tc_0_dt_0.01_num_1.dat');
V2=dlmread('V2_mean_5_sig_0_gNa_842.98_gCaT_10.1081_gCaS_12.5277_gA_112.336_gKCa_449.574_gK_281.021_gH_4.66558_gL_0.061134_gsyn_0_0_corr_0_tc_0_dt_0.01_num_1.dat');
figure;plot(0:0.1:5000,V1);hold on;plot(0:0.1:5000,V2)
xlim([4000 4500])

%2.
V1=dlmread('V1_mean_5_sig_0_gNa_859.992_gCaT_10.4957_gCaS_12.7299_gA_114.587_gKCa_458.594_gK_286.513_gH_4.58624_gL_0.12014_gsyn_0.573132_0.569253_corr_0_tc_0_dt_0.01_num_1.dat');
V2=dlmread('V2_mean_5_sig_0_gNa_853.375_gCaT_10.4_gCaS_12.5112_gA_113.798_gKCa_454.956_gK_284.416_gH_4.71101_gL_0.0813_gsyn_0.573132_0.569253_corr_0_tc_0_dt_0.01_num_1.dat');
figure;plot(0:0.1:5000,V1);hold on;plot(0:0.1:5000,V2)
xlim([4000 4500])
%without synapse
V1=dlmread('V1_mean_5_sig_0_gNa_859.992_gCaT_10.4957_gCaS_12.7299_gA_114.587_gKCa_458.594_gK_286.513_gH_4.58624_gL_0.12014_gsyn_0_0_corr_0_tc_0_dt_0.01_num_1.dat');
V2=dlmread('V2_mean_5_sig_0_gNa_853.375_gCaT_10.4_gCaS_12.5112_gA_113.798_gKCa_454.956_gK_284.416_gH_4.71101_gL_0.0813_gsyn_0_0_corr_0_tc_0_dt_0.01_num_1.dat');
figure;plot(0:0.1:5000,V1);hold on;plot(0:0.1:5000,V2)
xlim([4000 4500])
