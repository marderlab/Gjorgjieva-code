vec_test2 = [14973 12392 24635 3575 19668 43457 59045 10685 40404 18745 51845 35478 21858 51564 15902 65410 23438 53461 8111 45887 ...
             64973 50771 5820 20572 5884 14921 49320 61157 24631 34513 47279 15786 21760 51391 42901 41354 41379 47513 ... %spikers with different periods
             8226 10132 10155 10198 12699 12710 12717 15970 27317 27375 30165 33093 36348 41908 42299 ... %spikers with high g_H
             12570 12614 12622 17299 17316 17322 47830 51483 59185 59191 59196 59199 62986 62998 62999 63934 65091 65097 65866 ... %spikers with high g_A
             4993 5006 6386 6412 6426 45367 45377 45384 58368 58377 58384 59892 59899 62406 62419 63366 63375]; %spikers with low g_CaT and g_CaS and high g_KCa
ll=53; %1:6 8:9 11:14 16:21 24 26:29 33:34 36:45 47:53
vec_test2(ll)
format short g
tau = 100;
Esyn = -78;
Vhalf = -45;

vec = pairs(record2(vec_test2(ll)),:)
for ii=1:2
    g_Na(ii) = round(g(vec(ii),2)*10^6)/10^6;
    g_CaT(ii) = round(g(vec(ii),3)*10^6)/10^6;
    g_CaS(ii) = round(g(vec(ii),4)*10^6)/10^6;
    g_A(ii) = round(g(vec(ii),5)*10^6)/10^6;
    g_KCa(ii) = round(g(vec(ii),6)*10^6)/10^6;
    g_K(ii) = round(g(vec(ii),7)*10^6)/10^6;
    g_H(ii) = round(g(vec(ii),8)*10^6)/10^6;
end

gsyn_vec = [0.02:0.02:0.48];
vec_weird = zeros(length(gsyn_vec),length(gsyn_vec));
num_sp1 = nan(length(gsyn_vec),length(gsyn_vec));
num_sp2 = nan(length(gsyn_vec),length(gsyn_vec));
period1 = nan(length(gsyn_vec),length(gsyn_vec));
period2 = nan(length(gsyn_vec),length(gsyn_vec));
DC1 = nan(length(gsyn_vec),length(gsyn_vec));
DC2 = nan(length(gsyn_vec),length(gsyn_vec));
chi = nan(length(gsyn_vec),length(gsyn_vec));
corr_coef = nan(length(gsyn_vec),length(gsyn_vec));
phase_diff = nan(length(gsyn_vec),length(gsyn_vec));


bin_width = 10; %ms, for binning spike trains to compute correlations
offset = 0.5*10^4;
time_vec = [0:bin_width:2000*10^3];

num_sp1 = NaN*zeros(length(gsyn_vec),length(gsyn_vec));
num_sp2 = num_sp1;
DC1 = num_sp1;
DC2 = DC1;
period1 = DC1;
period2 = DC2;
chi = DC1;
corr_coef = DC1;
phase_diff = DC1;

vec_weird = zeros(length(gsyn_vec),length(gsyn_vec));
for g1=1:length(gsyn_vec)
    gsyn1 = gsyn_vec(g1);
    for g2=1:length(gsyn_vec)
        gsyn2 = gsyn_vec(g2);
        [g1 g2]
        fname = sprintf('sp_gNa_%g_%g_gCaT_%g_%g_gCaS_%g_%g_gA_%g_%g_gKCa_%g_%g_gK_%g_%g_gH_%g_%g_gL_0.01_0.01_gsyn_%g_%g_Vhalf_%g_tsyn_%g_Esyn_%g.dat',g_Na(1),g_Na(2),g_CaT(1),g_CaT(2),g_CaS(1),g_CaS(2),g_A(1),g_A(2),g_KCa(1),g_KCa(2),g_K(1),g_K(2),g_H(1),g_H(2),gsyn1,gsyn2,Vhalf,tau,Esyn);
        %cd data_half_centers/
        cd data_files_highH
        try
            sp = dlmread(fname);
        catch
            try
                structur = dir(sprintf('sp_gNa_%g_%g*gsyn_%g_%g*.dat',g_Na(1),g_Na(2),gsyn1,gsyn2));
                fname = structur.name;
                sp = dlmread(fname);
            catch
                try
                    structur = dir(sprintf('sp_*gCaS_%g_%g*gsyn_%g_%g*.dat',g_CaS(1),g_CaS(2),gsyn1,gsyn2));
                    fname = structur.name;
                    sp = dlmread(fname);
                catch
                    try 
                        structur = dir(sprintf('sp_*gA_%g_%g*gsyn_%g_%g*.dat',g_A(1),g_A(2),gsyn1,gsyn2));
                        fname = structur.name;
                        V = dlmread(fname);
                    catch
                        vec_weird(g1,g2) = 1; %if the whole spike file doesn't load
                    end
                end
            end
        end
        cd ..
        if (vec_weird(g1,g2)==0)
            sp1 = sp(sp(:,1)==0,2);
            sp2 = sp(sp(:,1)==1,2);
            if (length(sp1)>100)
                %bp1 = burst_params(sp1);
                %num_sp1(g1,g2) = bp1{1}(1);
                %eriod1(g1,g2) = bp1{1}(2);
                %DC1(g1,g2) = bp1{1}(3);
                [CV, num_sp, period_1, burst_dur, first_spike1, last_spike1] = get_train_properties_sp(sp1,0);
                num_sp1(g1,g2) = mean(num_sp);
                period1(g1,g2) = mean(period_1);
                DC1(g1,g2) = mean(burst_dur(1:end-1)./period_1);
            end
            if (length(sp2)>100) %sometimes there are very few spikes over 2000 sec
                %bp2 = burst_params(sp2);
                %num_sp2(g1,g2) = bp2{1}(1);
                %period2(g1,g2) = bp2{1}(2);
                %DC2(g1,g2) = bp2{1}(3);
                [CV, num_sp, period_2, burst_dur, first_spike2, last_spike2] = get_train_properties_sp(sp2,0);
                num_sp2(g1,g2) = mean(num_sp);
                period2(g1,g2) = mean(period_2);
                DC2(g1,g2) = mean(burst_dur(1:end-1)./period_2);
            end
            if (length(sp1)>100 & length(sp2)>100)
                if ((num_sp1(g1,g2)>1.01) & (num_sp2(g1,g2)>1.01))
                    disp('run chi');
                    ttt = get_exclusion_metric(sp1,sp2,5000);
                    chi(g1,g2) = -mean(ttt(~isnan(ttt)));
                end
                
                cell1_binned = histc(sp1-offset,time_vec); %
                cell2_binned = histc(sp2-offset,time_vec); %
                corr_coef(g1,g2) = xcov(cell1_binned,cell2_binned,0,'coef');
                
                %already computed first_spike and last_spike for both spikes trains
                if (length(sp1)<length(sp2)); first_spikee = first_spike1; first_spik = first_spike2; periodd = period_1;
                else first_spikee = first_spike2; first_spik = first_spike1; periodd = period_2; end
                
                phase_diff_vec = []; %must reset!
                for ts1=1:length(first_spikee)-1
                    [a,b] = min(abs(first_spikee(ts1)-first_spik));
                    if (a<periodd(ts1)) %if spike difference is much longer than period, don't count it, means gaps of no firing in data
                        pd = a/periodd(ts1);                  
                        phase_diff_vec = [phase_diff_vec pd];
                    end
                end
                phase_diff(g1,g2) = mean(phase_diff_vec);
            end
        end
    end
end
%sname = sprintf('props_pair_%d.mat',vec_test2(ll));
sname = sprintf('props_pair_highH_%d.mat',vec_test2(ll));
save(sname,'chi','corr_coef','num_sp1','num_sp2','period1','period2','DC1','DC2','phase_diff','vec_weird');

plot_vec = [gsyn_vec gsyn_vec(end)+0.01]-0.005;
generate_map;
len = length(gsyn_vec);

figure;subplot(3,2,1);pcolor(plot_vec,plot_vec,[[num_sp1; zeros(1,len)] zeros(len+1,1)]);caxis([min(min([num_sp1 num_sp2])) 1.001*max(max([num_sp1 num_sp2]))]);title('num spikes per burst');axis square;colorbar;set(gca,'ydir','normal');shading flat;set(gca,'color',0.5*[1 1 1]);xlim([min(plot_vec) max(plot_vec)]);ylim([min(plot_vec) max(plot_vec)])  
subplot(3,2,2);pcolor(plot_vec,plot_vec,[[num_sp2; zeros(1,len)] zeros(len+1,1)]);caxis([min(min([num_sp1 num_sp2])) 1.001*max(max([num_sp1 num_sp2]))]);title('num spikes per burst');axis square;colorbar;set(gca,'ydir','normal');shading flat;set(gca,'color',0.5*[1 1 1]);xlim([min(plot_vec) max(plot_vec)]);ylim([min(plot_vec) max(plot_vec)]) 
subplot(3,2,3);pcolor(plot_vec,plot_vec,[[period1; zeros(1,len)] zeros(len+1,1)]);caxis([min(min([period1 period2])) max(max([period1 period2]))]);title('period');axis square;colorbar;set(gca,'ydir','normal');shading flat;set(gca,'color',0.5*[1 1 1]);xlim([min(plot_vec) max(plot_vec)]);ylim([min(plot_vec) max(plot_vec)]) 
subplot(3,2,4);pcolor(plot_vec,plot_vec,[[period2; zeros(1,len)] zeros(len+1,1)]);caxis([min(min([period1 period2])) max(max([period1 period2]))]);title('period');axis square;colorbar;set(gca,'ydir','normal');shading flat;set(gca,'color',0.5*[1 1 1]);xlim([min(plot_vec) max(plot_vec)]);ylim([min(plot_vec) max(plot_vec)]) 
subplot(3,2,5);pcolor(plot_vec,plot_vec,[[DC1; zeros(1,len)] zeros(len+1,1)]);caxis([min(min([DC1 DC2])) max(max([DC1 DC2]))]);title('DC');axis square;colorbar;set(gca,'ydir','normal');shading flat;set(gca,'color',0.5*[1 1 1]);xlim([min(plot_vec) max(plot_vec)]);ylim([min(plot_vec) max(plot_vec)])  
subplot(3,2,6);pcolor(plot_vec,plot_vec,[[DC2; zeros(1,len)] zeros(len+1,1)]);caxis([min(min([DC1 DC2])) max(max([DC1 DC2]))]);title('DC');axis square;colorbar;set(gca,'ydir','normal');shading flat;set(gca,'color',0.5*[1 1 1]);xlim([min(plot_vec) max(plot_vec)]);ylim([min(plot_vec) max(plot_vec)])  
colormap(map);
set(gcf,'position',[0 0 600 700]);

figure;subplot(2,2,1);pcolor(plot_vec,plot_vec,[[chi; zeros(1,len)] zeros(len+1,1)]);axis square;colorbar;set(gca,'ydir','normal');title('burst exclusion metric');shading flat;set(gca,'color',0.5*[1 1 1]);xlim([min(plot_vec) max(plot_vec)]);ylim([min(plot_vec) max(plot_vec)]);%caxis([min(min(chi)) 1.001*max(max(chi))]); 
subplot(2,2,2);pcolor(plot_vec,plot_vec,[[corr_coef; zeros(1,len)] zeros(len+1,1)]);axis square;colorbar;set(gca,'ydir','normal');title('correlation coeffcient');shading flat;set(gca,'color',0.5*[1 1 1]); xlim([min(plot_vec) max(plot_vec)]);ylim([min(plot_vec) max(plot_vec)]) 
subplot(2,2,3);pcolor(plot_vec,plot_vec,[[DC1./DC2; zeros(1,len)] zeros(len+1,1)]);axis square;colorbar;set(gca,'ydir','normal');title('ratio of DC');shading flat; set(gca,'color',0.5*[1 1 1]);xlim([min(plot_vec) max(plot_vec)]);ylim([min(plot_vec) max(plot_vec)]) 
subplot(2,2,4);pcolor(plot_vec,plot_vec,[[phase_diff; zeros(1,len)] zeros(len+1,1)]);axis square;colorbar;set(gca,'ydir','normal');title('phase difference');shading flat;set(gca,'color',0.5*[1 1 1]);xlim([min(plot_vec) max(plot_vec)]);ylim([min(plot_vec) max(plot_vec)]) 
colormap(map);

%for some reason, before time 2000,000 seconds, this neuron is weird, but I
%don't have the voltage to prove it, look into it later. During the time
%for which voltage was recorded, the neuron is spiking at very high rate. 

%bp_temp = burst_params(sp(sp>2000000));
%bp(632,:)=bp_temp;