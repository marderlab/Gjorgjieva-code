function individual_PRCs(ll,ref_period1,ref_period2,time_burst_start1,time_burst_start2,gmax,Tdur1,Tdur2,TT,Vsyn)

load pairs_29Jan2016.mat
vec_test2 = [14973 12392 24635 3575 19668 43457 59045 10685 40404 18745 51845 35478 21858 51564 15902 65410 23438 53461 8111 45887 ...
             64973 50771 5820 20572 5884 14921 49320 61157 24631 34513 47279 15786 21760 51391 42901 41354 41379 47513 ... %spikers with different periods
             8226 10132 10155 10198 12699 12710 12717 15970 27317 27375 30165 33093 36348 41908 42299 ... %spikers with high g_H
             12570 12614 12622 17299 17316 17322 47830 51483 59185 59191 59196 59199 62986 62998 62999 63934 65091 65097 65866 ... %spikers with high g_A
             4993 5006 6386 6412 6426 45367 45377 45384 58368 58377 58384 59892 59899 62406 62419 63366 63375]; %spikers with low g_CaT and g_CaS and high g_KCa
    
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

tinit = 0;
tfinal = 3500; %ms
dt = 0.01;
tot_T = tinit:dt:tfinal;

ii = 1;
for frac=[0.01 0.02:0.02:0.98]%[0.01 0.1:0.1:0.9]
    Iapp = zeros(length(tot_T),2); %reset
    time1 = time_burst_start1 + frac*ref_period1;
    time_pulse = ceil(time1/dt):1:ceil((time1 + Tdur1)/dt);
    Iapp(time_pulse,1) = gmax; %time dependent vector, T x 2    
    time2 = time_burst_start2 + frac*ref_period2;
    time_pulse = ceil(time2/dt):1:ceil((time2 + Tdur2)/dt);
    Iapp(time_pulse,2) = gmax; %time dependent vector, T x 2    

    %figure;plot(tot_T,Iapp);xlim([2000 3000])
    [tot_T,V] = run_STG_model(g_Na,g_K,g_CaT,g_CaS,g_KCa,g_A,g_H,0,Iapp,Vsyn,0);
    V1 = V(:,1); V2 = V(:,2);
    Vth = 15; %mV threshold for detecting spikes
    dV1 = diff(V1); dV2 = diff(V2);
    x1 = find(dV1(1:end-1)>0 & dV1(2:end)<0);   %find where derivative switches sign     
    x2 = find(V1>Vth); 
    x3 = intersect(x1+1,x2);
    sp1 = tot_T(x3); %get spike times
    x1 = find(dV2(1:end-1)>0 & dV2(2:end)<0);   %find where derivative switches sign     
    x2 = find(V2>Vth); 
    x3 = intersect(x1+1,x2);
    sp2 = tot_T(x3); %get spike times

    [~, ~, period1, ~, first_spike1, ~] = get_train_properties_sp(sp1,0);
    [~, ~, period2, ~, first_spike2, ~] = get_train_properties_sp(sp2,0);

    time_burst_start_ind1 = find(first_spike1>2000,2); %get the bursts ignoring transient 2000 ms period
    time_burst_start_ind2 = find(first_spike2>2000,2); %get the bursts ignoring transient 2000 ms period
    PRC(ii,1) = (period1(time_burst_start_ind1(1)) - ref_period1)/ref_period1;
    PRC(ii,2) = (period2(time_burst_start_ind2(1)) - ref_period2)/ref_period2;
    PRC2(ii,1) = (period1(time_burst_start_ind1(2)) - ref_period1)/ref_period1;
    PRC2(ii,2) = (period2(time_burst_start_ind2(2)) - ref_period2)/ref_period2;
    ii = ii+1;
end
frac = [0.01 0.02:0.02:0.98];
namee = sprintf('PRCs_pair%g_gmax%g_Tdur%gper_Vsyn%g.mat',vec_test2(ll),gmax,TT*100,Vsyn);
save(namee,'PRC','PRC2','frac');