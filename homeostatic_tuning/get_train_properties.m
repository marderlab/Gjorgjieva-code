function [CV, num_sp, IBI, burst_dur, first_spike, last_spike] = get_train_properties(V1,t1,print_hist)

Vth = 15; %mV threshold for detecting spikes
%t1 = 0:0.1:5000; %voltage simulated with dt=0.01 but recorded at 0.1 resolution

dV1 = diff(V1);
x1 = find(dV1(1:end-1)>0 & dV1(2:end)<0);   %find where derivative switches sign     
x2 = find(V1>Vth); 
x3 = intersect(x1+1,x2);
sp1 = t1(x3); %get spike times

ISI = diff(sp1);  %get interspike interval   
CV = std(ISI)/mean(ISI);
vec = [0:0.1:500]; %bin for binning spike times
h1 = hist(ISI,vec); f1 = vec(h1>0);
if print_hist==1
    figure;plot(vec,h1)
    hold on;plot([sp_thr,sp_thr],[0 max(h1)],'k--')
end

sp_thr = 40;%ms - this is abitrary. Before I used mean(ISI), but it doesn't work for the second pair

n_bursts = 1;
if ~isempty(sp1) %make sure there's a nonzero number of spikes
    sp_thr = f1(end)/3.5;
    ii = 1; %take the first spike  
    n_bursts = 1; %counts the number of bursts
    num_sp(n_bursts) = 1; %first spike in that burst - note may not be a full burst
    first_spike(n_bursts) = sp1(ii);
    
    while length(sp1)>ii
        %for each burst
        if (sp1(ii+1)-sp1(ii) < sp_thr) %if spikes are close enough in time
            num_sp(n_bursts) = num_sp(n_bursts) + 1; %then add that spike to the current burst
        else
            last_spike(n_bursts) = sp1(ii); %prevous spike is the last spike in this burst
            n_bursts = n_bursts+1;
            num_sp(n_bursts) = 1; %new burst with 0 spikes
            first_spike(n_bursts) = sp1(ii+1); %next spike is first spike in the next burst            
        end
        ii = ii+1;
    end
    last_spike(n_bursts) = sp1(ii);
    %account for the fact that we may have not started counting at
    %the start of a burst...
    first_spike2 = first_spike(2:n_bursts-1);  
    last_spike2 = last_spike(2:n_bursts-1);
    num_sp = num_sp(2:n_bursts-1);
    IBI = diff(first_spike2)*0.001; %get into sec
    burst_dur = (last_spike2-first_spike2)*0.001; %sec
    burst_dur(burst_dur == 0) = 0.0001; %if the bursts consist of one spike each, i.e. tonic spiking
else
    num_sp = 0;
    IBI = 0;
    burst_dur = 0;
    first_spike = [];
    last_spike = [];
end

