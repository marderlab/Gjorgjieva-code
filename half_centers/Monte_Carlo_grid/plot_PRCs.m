vec_test2 = [14973 12392 24635 3575 19668 43457 59045 10685 40404 18745 51845 35478 21858 51564 15902 65410 23438 53461 8111 45887 ...
             64973 50771 5820 20572 5884 14921 49320 61157 24631 34513 47279 15786 21760 51391 42901 41354 41379 47513 ... %spikers with different periods
             8226 10132 10155 10198 12699 12710 12717 15970 27317 27375 30165 33093 36348 41908 42299 ... %spikers with high g_H
             12570 12614 12622 17299 17316 17322 47830 51483 59185 59191 59196 59199 62986 62998 62999 63934 65091 65097 65866 ... %spikers with high g_A
             4993 5006 6386 6412 6426 45367 45377 45384 58368 58377 58384 59892 59899 62406 62419 63366 63375]; %spikers with low g_CaT and g_CaS and high g_KCa
ll=2;

gmax_vec = [1 10 50 100 500]; %strength of synaptic conductance
Tdur_vec = [0.05 0.1 0.2 0.3 0.4]; %duration as percent of period
Vsyn_vec = [-78 0];  

for vv=1:length(Vsyn_vec)
    Vsyn = Vsyn_vec(vv); %for excitatory
    
    figure;
    for tt=1:length(Tdur_vec)        
        
        for gg=1:length(gmax_vec)
                gmax = gmax_vec(gg); %Tdur1 = 26; Tdur2 = Tdur1; %ms
                    
                namee = sprintf('PRCs_pair%g_gmax%g_Tdur%gper_Vsyn%g.mat',vec_test2(ll),gmax,Tdur_vec(tt)*100,Vsyn);
                load(namee);                
                
                subplot(2,5,tt); hold on; plot(frac,PRC(1:length(frac),1),'o-');ylim([-1 1])
                subplot(2,5,tt+5); hold on; plot(frac,PRC(1:length(frac),2),'o-');ylim([-1 1])
        end
    end
end
                
%%


for vv=1:length(Vsyn_vec)
    Vsyn = Vsyn_vec(vv); %for excitatory
    
    figure;
    for tt=1:length(Tdur_vec)        
        
        for gg=1:length(gmax_vec)
                gmax = gmax_vec(gg); %Tdur1 = 26; Tdur2 = Tdur1; %ms
                    
                namee = sprintf('PRCs_HCO_pair%g_gmax%g_Tdur%gper_Vsyn%g.mat',vec_test2(ll),gmax,Tdur_vec(tt)*100,Vsyn);
                try load(namee);                
                
                subplot(2,5,tt); hold on; plot(frac,PRC_n2(1:length(frac),1),'o-');ylim([-1 1])
                subplot(2,5,tt+5); hold on; plot(frac,PRC_n2(1:length(frac),2),'o-');ylim([-1 1])
                catch;end
        end
    end
end
      
                
                
                