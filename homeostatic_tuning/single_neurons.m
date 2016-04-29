gA_vec = 0:10:100;
I_vec = 0:0.5:10;
g_syn = 0.3;

Vth = 15;
t1 = 0:0.1:5000; %the total time over which voltage was recorded
for aa=1:length(gA_vec)
    gA = gA_vec(aa);
    for vv = 1:length(I_vec)
        I = I_vec(vv);
        
        V1=dlmread(sprintf('V2_mean_%g_sig_0_gNa_227.052_gCaT_2.7781_gCaS_3.3469_gA_%g_gKCa_121.124_gK_75.649_gH_1.3194_gL_0.1631_gsyn_%g_%g_corr_0_tc_0_dt_0.01_num_1.dat',I,gA,g_syn,g_syn));
        %V2=dlmread(sprintf('V2_mean_%g_sig_0_gNa_227.052_gCaT_2.7781_gCaS_3.3469_gA_%g_gKCa_121.124_gK_75.649_gH_1.3194_gL_0.1631_gsyn_%g_%g_corr_0_tc_0_dt_0.01_num_0.dat',I,gA,g_syn,g_syn));
        %figure;plot(t1,V1);
        %hold on;plot(t1,V2);
        %xlim([4500 5000]); title([gA,I]);pause(1);close
        dV1 = diff(V1);
        x1 = find(dV1(1:end-1)>0 & dV1(2:end)<0);   %find where derivative switches sign     
        x2 = find(V1>Vth); 
        x3 = intersect(x1+1,x2);
        %figure;plot(t1,V1);hold on;plot(t1(x3),15,'ro')
        sp1 = t1(x3);
        %figure;hist(diff(sp1),20);
        
        [CV, num_sp, IBI, burst_dur, first_spike, last_spike] = get_train_properties(V1,t1,0);
        %pause(1);close
        try
            num_spikes(aa,vv) = mean(num_sp);
        catch
            num_spikes(aa,vv) = 0;
        end
        try
            burst_freq(aa,vv) = 1/mean(IBI); %Hz
        catch
            burst_freq(aa,vv) = 1/0; %Hzk
        end
        %sp_freq(aa,vv) = mean(num_sp./burst_dur); %Hz
        %DC(aa,vv) = mean(burst_dur(1:end-1)./IBI); %duty cycle
    end
end
%max_b = max(max(burst_freq(burst_freq<Inf)));
%min_b = min(min(burst_freq));

burst_freq(burst_freq==Inf) = -1;%max_b + (max_b-min_b)/63;

figure;imagesc(I_vec,gA_vec,burst_freq);set(gcf,'Position', [100, 100, 600, 280]);colorbar;set(gca,'ydir','normal')
figure;imagesc(I_vec,gA_vec,num_spikes);set(gcf,'Position', [100, 100, 600, 280]);colorbar;set(gca,'ydir','normal')