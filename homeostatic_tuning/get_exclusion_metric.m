function chi = get_exclusion_metric(sp1,sp2,Tp)
%get burst exclusion metric for two spike trains
%Tp indicates over what period is the metrix computed, as a running average

%remove first 5000 sec, transient period
%normalize so that spike trains start at 0
spp1 = sp1(sp1>5000)-5000; 
spp2 = sp2(sp2>5000)-5000;

%over what time period?
%Tp = 10000; %10 sec
Tslide = 100; %ms
Ttot = 2000000; %2000 sec

%compute the burst exclusion metric for multiple shorter intervals of
%length 5000 which move with a sliding window of 100
count = 0;
for ii=Tp:Tslide:Ttot-Tp
    %normalize each time so that the new spike trains starts at 0 and ends
    %at Tp (total duration of each short spike train interval
    sp11 = spp1(spp1>ii & spp1<ii+Tp)-ii; 
    sp12 = spp2(spp2>ii & spp2<ii+Tp)-ii;
    count = count+1;
    chi2(count) = burst_exclusion(sp11,sp12,Tp);
end

chi = mean(chi2);
%figure;plot(Tp:Tslide:Ttot-Tp,chi2,'o-')
