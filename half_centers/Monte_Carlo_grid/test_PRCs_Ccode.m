vec_test2 = [14973 12392 24635 3575 19668 43457 59045 10685 40404 18745 51845 35478 21858 51564 15902 65410 23438 53461 8111 45887 ...
             64973 50771 5820 20572 5884 14921 49320 61157 24631 34513 47279 15786 21760 51391 42901 41354 41379 47513 ... %spikers with different periods
             8226 10132 10155 10198 12699 12710 12717 15970 27317 27375 30165 33093 36348 41908 42299 ... %spikers with high g_H
             12570 12614 12622 17299 17316 17322 47830 51483 59185 59191 59196 59199 62986 62998 62999 63934 65091 65097 65866 ... %spikers with high g_A
             4993 5006 6386 6412 6426 45367 45377 45384 58368 58377 58384 59892 59899 62406 62419 63366 63375]; %spikers with low g_CaT and g_CaS and high g_KCa
ll=1; %all ll are ok for the first data set load pairs_29Jan2016.mat
%1:6 8:9 11:14 16:21 24 26:29 33:34 36:45 47:53 %for the second data set load pairs_Feb16.mat
load pairs_29Jan2016.mat
format short g
tau = 100;
Esyn = -78;
Vhalf = -45;

tinit = 0;
tfinal = 5000; %ms
dt = 0.1;
tot_T = tinit:dt:tfinal;
 

ll = 10; %1:53 %go through all the pairs
vec = pairs(record2(vec_test2(ll)),:);
for ii=1:2
    g_Na(ii) = round(g(vec(ii),2)*10^6)/10^6;
    g_CaT(ii) = round(g(vec(ii),3)*10^6)/10^6;
    g_CaS(ii) = round(g(vec(ii),4)*10^6)/10^6;
    g_A(ii) = round(g(vec(ii),5)*10^6)/10^6;
    g_KCa(ii) = round(g(vec(ii),6)*10^6)/10^6;
    g_K(ii) = round(g(vec(ii),7)*10^6)/10^6;
    g_H(ii) = round(g(vec(ii),8)*10^6)/10^6;
end

%first run it without input
V_18745_0 = dlmread('V_in1_312_in2_1125_gsyn_0_0.dat');
V1 = V_18745_0(V_18745_0(:,1)==0,2); 
V2 = V_18745_0(V_18745_0(:,1)==1,2);
Vth = 15; %mV threshold for detecting spikes
dV1 = diff(V1); 
dV2 = diff(V2);
x1 = find(dV1(1:end-1)>0 & dV1(2:end)<0);   %find where derivative switches sign     
x2 = find(V1>Vth); 
x3 = intersect(x1+1,x2);
sp1 = tot_T(x3); %get spike times
x1 = find(dV2(1:end-1)>0 & dV2(2:end)<0);   %find where derivative switches sign     
x2 = find(V2>Vth); 
x3 = intersect(x1+1,x2);
sp2 = tot_T(x3); %get spike times
[CV1, num_sp1, period1, burst_dur1, first_spike1, last_spike1] = get_train_properties_sp(sp1,0);
[CV2, num_sp2, period2, burst_dur2, first_spike2, last_spike2] = get_train_properties_sp(sp2,0);

%disp('how many bursts?')
%[length(first_spike1) length(first_spike2)]
time_burst_start1 = first_spike1(find(first_spike1>2000,1)) %get the bursts ignoring transient 2000 ms period
time_burst_start2 = first_spike2(find(first_spike2>2000,1)) %get the bursts ignoring transient 2000 ms period
ref_period1 = period1(find(first_spike1>tfinal*2/3,1))
ref_period2 = period2(find(first_spike2>tfinal*2/3,1))

%now run it with input
V_18745 = dlmread('V_in1_312_in2_1125_gsyn_0.07_0.09.dat');
V1 = V_18745(V_18745(:,1)==0,2); 
V2 = V_18745(V_18745(:,1)==1,2);
Vth = 15; %mV threshold for detecting spikes
dV1 = diff(V1); 
dV2 = diff(V2);
x1 = find(dV1(1:end-1)>0 & dV1(2:end)<0);   %find where derivative switches sign     
x2 = find(V1>Vth); 
x3 = intersect(x1+1,x2);
sp1 = tot_T(x3); %get spike times
x1 = find(dV2(1:end-1)>0 & dV2(2:end)<0);   %find where derivative switches sign     
x2 = find(V2>Vth); 
x3 = intersect(x1+1,x2);
sp2 = tot_T(x3); %get spike times
[CV1, num_sp1, period1, burst_dur1, first_spike1, last_spike1] = get_train_properties_sp(sp1,0);
[CV2, num_sp2, period2, burst_dur2, first_spike2, last_spike2] = get_train_properties_sp(sp2,0);

%disp('how many bursts?')
%[length(first_spike1) length(first_spike2)]
time_burst_start1 = first_spike1(find(first_spike1>2000,1)) %get the bursts ignoring transient 2000 ms period
time_burst_start2 = first_spike2(find(first_spike2>2000,1)) %get the bursts ignoring transient 2000 ms period
ref_period1 = period1(find(first_spike1>tfinal*2/3,1))
ref_period2 = period2(find(first_spike2>tfinal*2/3,1))