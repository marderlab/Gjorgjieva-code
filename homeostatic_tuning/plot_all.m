Vth = 15; %mV threshold for detecting spikes
t1 = 0:0.1:5000; %voltage simulated with dt=0.01 but recorded at 0.1 resolution
gA_vec = 0:10:100;
I_vec = 0:0.5:10;
g_syn = 0;

for aa=1:length(gA_vec)
    gA = gA_vec(aa);
    for vv = 1:length(I_vec)
        I = I_vec(vv);
        clear first_spike1 last_spike1 num_sp1 burst_dur IBI

        V1=dlmread(sprintf('V1_mean_%g_sig_0_gNa_227.052_gCaT_2.7781_gCaS_3.3469_gA_%g_gKCa_121.124_gK_75.649_gH_1.3194_gL_0.1631_gsyn_%g_%g_corr_0_tc_0_dt_0.01_num_0.dat',I,gA,g_syn,g_syn));
        %V2=dlmread(sprintf('V2_mean_%g_sig_0_gNa_227.052_gCaT_2.7781_gCaS_3.3469_gA_%g_gKCa_121.124_gK_75.649_gH_1.3194_gL_0.1631_gsyn_%g_%g_corr_0_tc_0_dt_0.01_num_0.dat',I,gA,g_syn,g_syn));
        %figure;plot(0:0.1:5000,V1);
        %hold on;plot(0:0.1:5000,V2);
        %xlim([4500 5000]); title([gA,I]);pause(1);closex
        dV1 = diff(V1);
        x1 = find(dV1(1:end-1)>0 & dV1(2:end)<0);   %find where derivative switches sign     
        x2 = find(V1>Vth); 
        x3 = intersect(x1+1,x2);
        %figure;plot(t1,V1);hold on;plot(t1(x3),15,'ro')
        sp1 = t1(x3);
        %figure;hist(diff(sp1),20);
        
        ISI1 = diff(sp1);       
        vec = [0:0.1:500];
        h1 = hist(ISI1,vec); f1 = vec(h1>0);
        n_bursts = 1;
        if ~isempty(sp1) %make sure there's a nonzero number of spikes
            ii = 1; %take the first spike  
            n_bursts = 1; %counts the number of bursts
            num_sp1(n_bursts) = 1; %first spike in that burst - note may not be a full burst
            
            while length(sp1)>ii
                %for each burst
                if (sp1(ii+1)-sp1(ii) < f1(end)/3) %if spikes are close enough in time
                    num_sp1(n_bursts) = num_sp1(n_bursts) + 1; %then add that spike to the current burst
                else
                    first_spike1(n_bursts) = sp1(ii+1); %next spike is first spike in the next burst
                    last_spike1(n_bursts) = sp1(ii); %prevous spike is the last spike in this burst
                    n_bursts = n_bursts+1;
                    num_sp1(n_bursts) = 1; %new burst with 0 spikes
                end
                ii = ii+1;
            end
            %account for the fact that we may have not started counting at
            %the start of a burst...
            first_spike1 = first_spike1(1:n_bursts-2);  
            last_spike1 = last_spike1(2:n_bursts-1);
            num_sp1 = num_sp1(1:n_bursts-2);
            num_spikes1(aa,vv) = mean(num_sp1);
            IBI = diff(first_spike1)*0.001; %get into sec
            burst_freq1(aa,vv) = 1/mean(IBI); %Hz
            burst_dur = (last_spike1-first_spike1)*0.001; %sec
            burst_dur(burst_dur == 0) = 0.0001;
            sp_freq1(aa,vv) = mean(num_sp1./burst_dur); %Hz
            DC1(aa,vv) = mean(burst_dur(1:end-1)./IBI); %duty cycle
        else
            burst_freq1(aa,vv) = 0;
            sp_freq1(aa,vv) = 0;
            DC1(aa,vv) = 0;
        end
    end
end