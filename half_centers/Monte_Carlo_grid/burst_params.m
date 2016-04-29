% bp = [num_sp_burst, period, DC]
function bp = burst_params(st)
%st = spike times in ms
%based on code by Tim O'Leary from Neuron 2014 paper

ISI = diff(st); 
mid = (max(ISI) + min(ISI))/2;
fastISI = mean(ISI(ISI < mid)); %ISI
freq = 1000/fastISI;
IBI = (ISI(ISI >= mid)); %mean(IBI)
first_spikes = ([0; diff(ISI < mid)] > 0);
first_spike_times = st(first_spikes);
last_spikes = (diff(ISI < mid) > 0);
last_spike_times = st(last_spikes);
first_spike_times = first_spike_times(1:end-1);
last_spike_times = last_spike_times(2:end);

% num_sp = [];
% for ii=1:length(first_spike_times)
%     num_sp = [num_sp length(st(st>=first_spike_times(ii)-eps & st<=last_spike_times(ii)+eps))];
% end
%a more efficient way to implement this is:
%total # sp/total # bursts
%but... if there's a pattern that's 2-3-2-3-2... then it won't give us the
%vector of all num_sp in each burst, just the mean
num_sp = length(st)/length(first_spike_times+1);

%spiking neuron
if (length(first_spike_times)<5 | abs(fastISI-mean(IBI))<1)
    %error('not enough bursts (< 5)');
    bp{1} = [1 fastISI 0];
    bp{2} = st;
    return;
end

%figure;plot(num_sp,'o')
num_sp_burst = mean(num_sp);

period = mean(diff(first_spike_times)); %burst period
DC = (period - mean(IBI))/period; %duty cycle
fastISI = mean(ISI(ISI < mid)); %ISI

% everything up to here in samples (not time units)

% %%%%%%%%%% find threshold
% Commented by JG
% seglength = floor(mean(IBI));
% thrvals = zeros(length(first_spike_times-2)-1,1);
% mins = thrvals;
% for i=2:length(first_spike_times)
%     try
%         vseg = v(first_spike_times(i)-seglength:first_spike_times(i));
% 
%         mins(i-1) = min(vseg);
% 
%         k = deriv(vseg,1e-6,2)./((1+deriv(vseg,1e-6).^2).^1.5);
% 
%         %plot(vseg);
%         %hold on;
%         %plot(k);
% 
%         thridx = findhump(k,'backwards','');
%         thrvals(i-1) = vseg(thridx);
%     catch
%         'error'
%         bp = nan(1,5);
%         return;
%     end
% end

%%%%%%%%%%% find spike heights
% Commented by JG
% spheight = zeros(length(st)-1,1);
% 
% for i=1:length(st)-2
%     %i
%     %length(st)
%     vseg = v(st(i):st(i) + floor(fastISI));
%     maxidx = findhump(vseg,'forwards','');
%     if ~isempty(maxidx)
%         spheight(i) = vseg(maxidx);
%     else
%         spheight(i) = NaN;
%     end
% end

%plot(v,'b');
%hold on;
%plot(st(1:end-1),spheight,'*k');
%plot(first_spike_times,mean(spheight),'*r');

% Commented by JG
% slomax = mean(thrvals);
% sh = mean(spheight) - slomax;
% sloamp = slomax - mean(mins);

%plot(v,'b');
%v(v>slomax) = slomax;
%hold on;
%plot(v,'r');

bp{1} = [num_sp_burst period DC];
bp{2} = first_spike_times;


% Commented by JG
% function s = spikeTimes(v,t)
% 
% s = find(diff(v >= t) > 0);
% end


