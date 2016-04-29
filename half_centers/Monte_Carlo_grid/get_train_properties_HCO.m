function [CV, num_sp, period, burst_dur, first_spike2, last_spike2, sp_thr] = get_train_properties_HCO(sp1,print_hist,sp_thr)

%provided spike times here, not just voltage!

ISI = diff(sp1);  %get interspike interval   
CV = std(ISI)/mean(ISI);
vec = [0:0.1:500]; %bin for binning spike times
h1 = hist(ISI,vec); f1 = vec(h1>0);

%sp_thr = 40;%ms - this is abitrary. Before I used mean(ISI), but it doesn't work for the second pair
ISI = diff(sp1); 
%sorted_ISI = sort(ISI);

if nargin==2
   sp_thr = (quantile(ISI,0.9) + min(ISI))/2;
    %sp_thr = (max(ISI)*1.1 + min(ISI))/2;
    if (abs(sp_thr-max(ISI))<10 | abs(sp_thr-min(ISI))<10)
        sp_thr = 0.99*min(ISI);
    end
end
if print_hist==1
    figure;hist(ISI,20);%bar(vec,h1)
    hold on;plot([sp_thr,sp_thr],[0 max(h1)],'k--')
    quantile(ISI,0.95)
    sp_thr
end
n_bursts = 1;
if ~isempty(sp1) %make sure there's a nonzero number of spikes
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
            n_bursts = n_bursts+1; %increase number of bursts
            num_sp(n_bursts) = 1; %new burst with 1 spike
            first_spike(n_bursts) = sp1(ii+1); %next spike is first spike in the next burst            
        end
        ii = ii+1;
    end
    last_spike(n_bursts) = sp1(ii);
    %account for the fact that we may have not started counting at
    %the start of a burst OR edned at the end of a burst
    first_spike2 = first_spike(2:n_bursts-1);  
    last_spike2 = last_spike(2:n_bursts-1);
    num_sp = num_sp(2:n_bursts-1);
    period = diff(first_spike2); %ms
    burst_dur = last_spike2-first_spike2; %ms
    burst_dur(burst_dur == 0) = 0.0001; %if the bursts consist of one spike each, i.e. tonic spiking
else
    num_sp = NaN;
    period = NaN;
    burst_dur = NaN;
    first_spike2 = [];
    last_spike2 = [];
end

