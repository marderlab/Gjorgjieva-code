clear num_spikes1 num_spikes2 burst_freq1 burst_freq2 sp_freq1 sp_freq2 DC1 DC2
mu_vec = [0:0.1:10];
g_vec = [0:0.001:0.02];

vec = [0:0.1:500]; %for binning of the ISI
n_bursts = 50; %number of bursts to use in the analysis, determine oscillation properties

for mm=1:length(mu_vec)
    mu = mu_vec(mm);
    for gg=1:length(g_vec)
        g_syn = g_vec(gg);
   
        [mu g_syn]
        try
            sp1 = dlmread(sprintf('dat_files/spikes1_mean_%g_sig_0_gNa_227.052_gCaT_2.7781_gCaS_3.3469_gA_30.4321_gKCa_121.124_gK_75.649_gH_1.3194_gL_0.1631_gsyn_%g_%g_corr_0_tc_1_dt_0.01_num_2.dat',mu,g_syn,g_syn));
            %sp1 = dlmread(sprintf('two_diff_pair/spikes1_mean_%g_sig_0_gNa_227.052_gCaT_2.7781_gCaS_3.3469_gA_30.4321_gKCa_121.124_gK_75.649_gH_1.3194_gL_0.1631_gsyn_%g_%g_corr_0_tc_1_dt_0.01_num_2.dat',mu,g_syn,g_syn));            
        catch
            sp1 = [];
        end
        try 
            sp2 = dlmread(sprintf('dat_files/spikes2_mean_%g_sig_0_gNa_227.052_gCaT_2.7781_gCaS_3.3469_gA_30.4321_gKCa_121.124_gK_75.649_gH_1.3194_gL_0.1631_gsyn_%g_%g_corr_0_tc_1_dt_0.01_num_2.dat',mu,g_syn,g_syn));
            %sp2 = dlmread(sprintf('two_diff_pair/spikes2_mean_%g_sig_0_gNa_158.687_gCaT_2.0942_gCaS_2.3907_gA_21.3218_gKCa_84.7673_gK_52.9122_gH_0.8532_gL_0.0561_gsyn_%g_%g_corr_0_tc_1_dt_0.01_num_2.dat',mu,g_syn,g_syn));
        catch
            sp2 = [];
        end
        ISI1 = diff(sp1); ISI2 = diff(sp2);
        h1 = hist(ISI1,[0:0.1:500]); f1 = vec(h1>0);
        h2 = hist(ISI2,[0:0.1:500]); f2 = vec(h2>0);

%         d1 = diff(f1);
%         count = 0;
%         ii = 1;
%         while (ii<length(d1)+1)
%             if (d1(ii)==1)
%                 if (ii==1)
%                     count = count+1;
%                     f1_new(count) = f1(ii);
%                 end
%                 if (ii<length(d1) && d1(ii+1)==1)
%                     count = count+1;
%                     f1_new(count) = f1(ii+1);
%                 end
%                 if (ii==length(d1))
%                     count = count+1;
%                     f1_new(count) = f1(ii+1);
%                 end
%             end
%             if (d1(ii)==0)
%                 numm = 1;
%                 summ = f1(ii);
%                 while (ii<length(d1)+1 && d1(ii)==0)
%                     numm = numm+1;
%                     summ = summ + f1(ii+1);
%                     ii = ii+1;
%                 end
%                 ii = ii-1;
%                 count = count+1;
%                 f1_new(count) = summ/numm;
%             end
%             ii = ii+1;
%         end
        %make sure we consider spikes after a stationary state has been reached
        xx = find(sp1>40000);
        if length(xx)>0 %make sure there's a nonzero number of spikes
            ii = xx(1); %take the first spike  
            for count = 1:n_bursts+1 %look at a given number of bursts (n_bursts)...  
                %for each burst
                num_sp1(count) = 1; %first spike in that burst - note may not be a full burst
                while (sp1(ii+1)-sp1(ii) < f1(end)/3) %if spikes are close enough in time
                    ii = ii+1;
                    num_sp1(count) = num_sp1(count) + 1; %then add that spike to the current burst
                end
                first_spike1(count) = sp1(ii+1); %next spike is first spike in the next burst
                last_spike1(count) = sp1(ii); %prevous spike is the last spike in this burst
                ii = ii+1;
            end
            %account for the fact that we may have not started counting at
            %the start of a burst...
            first_spike1 = first_spike1(1:n_bursts);
            last_spike1 = last_spike1(2:n_bursts+1);
            num_sp1 = num_sp1(2:end);
            num_spikes1(mm,gg) = mean(num_sp1);
            IBI = diff(first_spike1)*0.001; %get into sec
            burst_freq1(mm,gg) = 1/mean(IBI); %Hz
            burst_dur = (last_spike1-first_spike1)*0.001; %sec
            burst_dur(burst_dur == 0) = 0.0001;            
            sp_freq1(mm,gg) = mean(num_sp1./burst_dur); %Hz
            DC1(mm,gg) = mean(burst_dur(1:end-1)./IBI); %duty cycle
        else
            num_spikes1(mm,gg) = 0;
            burst_freq1(mm,gg) = 0;
            sp_freq1(mm,gg) = 0;
            DC1(mm,gg) = 0;
            phi(mm,gg) = 0; %initialize
        end
        
        xx = find(sp2>40000);
        if length(xx)>0
            ii = xx(1);      
            for count = 1:n_bursts+1
                num_sp2(count) = 1;
                while (sp2(ii+1)-sp2(ii) < f2(end)/3)
                    ii = ii+1;
                    num_sp2(count) = num_sp2(count) + 1;
                end
                first_spike2(count) = sp2(ii+1);
                last_spike2(count) = sp2(ii);
                ii = ii+1;
            end
            first_spike2 = first_spike2(1:n_bursts);
            last_spike2 = last_spike2(2:n_bursts+1);
            num_sp2 = num_sp2(2:end);
            num_spikes2(mm,gg) = mean(num_sp2);
            IBI = diff(first_spike2)*0.001; %sec
            burst_freq2(mm,gg) = 1/mean(IBI);
            burst_dur = (last_spike2-first_spike2)*0.001;
            burst_dur(burst_dur == 0) = 0.0001;
            sp_freq2(mm,gg) = mean(num_sp2./burst_dur);
            DC2(mm,gg) = mean(burst_dur(1:end-1)./IBI);
        else
            num_spikes2(mm,gg) = 0;
            burst_freq2(mm,gg) = 0;
            sp_freq2(mm,gg) = 0;
            DC2(mm,gg) = 0;
            phi(mm,gg) = 0;
        end
        
        phi(mm,gg) = min(abs(mean(first_spike1 - first_spike2)*0.001/mean(IBI)),1-abs(mean(first_spike1 - first_spike2)*0.001/mean(IBI)));
            
    end
end
