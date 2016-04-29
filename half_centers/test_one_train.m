clear num_spikes1 num_spikes2 burst_freq1 burst_freq2 sp_freq1 sp_freq2 DC1 DC2

mm = 1;
gg = 1;
mu = 0.9;
g_syn = 0.01;
sp1 = dlmread(sprintf('dat_files/spikes1_mean_%g_sig_0_gNa_227.052_gCaT_2.7781_gCaS_3.3469_gA_30.4321_gKCa_121.124_gK_75.649_gH_1.3194_gL_0.1631_gsyn_%g_%g_corr_0_tc_1_dt_0.01_num_2.dat',mu,g_syn,g_syn));

ISI1 = diff(sp1); 
h1 = hist(ISI1,[0:0.1:500]); f1 = vec(h1>0);
n_bursts = 30;
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
    burst_freq1(mm,gg) = 0;
    sp_freq1(mm,gg) = 0;
    DC1(mm,gg) = 0;
end