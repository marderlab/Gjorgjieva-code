pairs = nchoosek(1:2000,2);

num_pairs = size(pairs,1);

%%
%record with th_num_sp=0.5;th_period=5; th_DC=0.05; no thr on V_hyp -> 36308 pairs
%record_V with th_num_sp=0.5;th_period=5; th_DC=0.05;th_V=0.5 mV -> 4922 pairs
%record2 with th_num_sp=0.99;th_period=20; th_DC=0.1;th_V=5 mV -> 138576 pairs
%record2 with th_num_sp=0.99;th_period=20; th_DC=0.1;th_V=5 mV; MSE(g1,g2)>median -> 66211 pairs

record2 = [];
th_num_sp = 0.99; %similar number of spikes per burst
th_period = 20; %similar perid within this many ms
th_DC = 0.1; %similar duty cycle
th_V = 5; %similar hyperpol V
med_MSE = median(MSE);
mean_MSE = mean(MSE);
med_L1 = median(L1);
mean_L1 = mean(L1);
for ii=1:num_pairs
    diff_bp1 = abs(bp(pairs(ii,1),1)-bp(pairs(ii,2),1));
    diff_bp2 = abs(bp(pairs(ii,1),2)-bp(pairs(ii,2),2));
    diff_bp3 = abs(bp(pairs(ii,1),3)-bp(pairs(ii,2),3));
    diff_bp4 = abs(most_hyper_V(pairs(ii,1)) - most_hyper_V(pairs(ii,2)));
    diff_g = mean((g(pairs(ii,1),2:8)-g(pairs(ii,2),2:8)).^2);
    if (diff_bp1<th_num_sp & diff_bp2<th_period & diff_bp3<th_DC & diff_bp4<th_V & diff_g>med_MSE)
        record2 = [record2 ii];
    end
end

%%
figure;subplot(2,2,1);plot(bp(pairs(record2),1),'bo-');xlim([0 length(record2)]);title('num spikes/burst');
subplot(2,2,2);plot(bp(pairs(record2),2),'ro-');xlim([0 length(record2)]);title('burst period (ms)')
subplot(2,2,3);plot(bp(pairs(record2),3),'ko-');xlim([0 length(record2)]);title('DC');
subplot(2,2,4);plot(most_hyper_V(pairs(record2)),'o-','color',[0.5 0.5 0.5]);xlim([0 length(record2)]);title('V_{hyper} (mV)');
%%
%to make sure that not all the similar pairs that have been chosen are
%spiking neurons...
%figure;plot(bp(pairs(record),1))

for num_rec = 26453%1:length(record)

ii = pairs(record(num_rec),1);
g_Na = round(g(ii,2)*10^6)/10^6;
g_CaT = round(g(ii,3)*10^6)/10^6;
g_CaS = round(g(ii,4)*10^6)/10^6;
g_A = round(g(ii,5)*10^6)/10^6;
g_KCa = round(g(ii,6)*10^6)/10^6;
g_K = round(g(ii,7)*10^6)/10^6;
g_H = round(g(ii,8)*10^6)/10^6;
g_L = 0.01;
fname = sprintf('sp_gNa_%g_gCaT_%g_gCaS_%g_gA_%g_gKCa_%g_gK_%g_gH_%g_gL_%g.dat',g_Na,g_CaT,g_CaS,g_A,g_KCa,g_K,g_H,g_L);
fname2 = sprintf('V_gNa_%g_gCaT_%g_gCaS_%g_gA_%g_gKCa_%g_gK_%g_gH_%g_gL_%g.dat',g_Na,g_CaT,g_CaS,g_A,g_KCa,g_K,g_H,g_L);

cd data_files_from_cluster/
try
    sp = dlmread(fname);
    V = dlmread(fname2);
catch
    try
        structur = dir(sprintf('sp_gNa_%g_*.dat',g_Na));
        fname = structur.name;
        sp = dlmread(fname);
        structur = dir(sprintf('V_gNa_%g_*.dat',g_Na));
        fname2 = structur.name;
        V = dlmread(fname2);
    catch
        try
            structur = dir(sprintf('sp_*gCaS_%g_*.dat',g_CaS));
            fname = structur.name;
            sp = dlmread(fname);
            structur = dir(sprintf('V_*gCaS_%g_*.dat',g_CaS));
            fname2 = structur.name;
            V = dlmread(fname2);
        catch
            vec_weird = [vec_weird ii];
        end
    end
end
cd ..
sp1 = sp;
V1 = V;

%%%%%% second pair
ii = pairs(record(num_rec),2);
g_Na = round(g(ii,2)*10^6)/10^6;
g_CaT = round(g(ii,3)*10^6)/10^6;
g_CaS = round(g(ii,4)*10^6)/10^6;
g_A = round(g(ii,5)*10^6)/10^6;
g_KCa = round(g(ii,6)*10^6)/10^6;
g_K = round(g(ii,7)*10^6)/10^6;
g_H = round(g(ii,8)*10^6)/10^6;
g_L = 0.01;
fname = sprintf('sp_gNa_%g_gCaT_%g_gCaS_%g_gA_%g_gKCa_%g_gK_%g_gH_%g_gL_%g.dat',g_Na,g_CaT,g_CaS,g_A,g_KCa,g_K,g_H,g_L);
fname2 = sprintf('V_gNa_%g_gCaT_%g_gCaS_%g_gA_%g_gKCa_%g_gK_%g_gH_%g_gL_%g.dat',g_Na,g_CaT,g_CaS,g_A,g_KCa,g_K,g_H,g_L);

cd data_files_from_cluster/
try
    sp = dlmread(fname);
    V = dlmread(fname2);
catch
    try
        structur = dir(sprintf('sp_gNa_%g_*.dat',g_Na));
        fname = structur.name;
        sp = dlmread(fname);
        structur = dir(sprintf('V_gNa_%g_*.dat',g_Na));
        fname2 = structur.name;
        V = dlmread(fname2);
    catch
        try
            structur = dir(sprintf('sp_*gCaS_%g_*.dat',g_CaS));
            fname = structur.name;
            sp = dlmread(fname);
            structur = dir(sprintf('V_*gCaS_%g_*.dat',g_CaS));
            fname2 = structur.name;
            V = dlmread(fname2);
        catch
            vec_weird = [vec_weird ii];
        end
    end
end
cd ..
sp2 = sp;
V2 = V;

t = 0:0.1:5000;
figure(100);plot(t,V1); hold on;plot(t,V2,'r--');xlim([3000 5000])
title(sprintf('%g of %g',num_rec,length(record)));
pause(0.2); hold off
end