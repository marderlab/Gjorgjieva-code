ygNa_min = 800; gNa_max = 1200;
gCaT_min = 0; gCaT_max = 6;
gCaS_min = 0; gCaS_max = 12;
gA_min = 20; gA_max = 130;
gKCa_min = 20; gKCa_max = 140;
gK_min = 90; gK_max = 120;
gH_min = 1.5; gH_max = 3;

Num_sims = 200;
for ii = 3:Num_sims
    ii
    g_Na(ii) = rand*(gNa_max-gNa_min) + gNa_min;
    g_CaT(ii) = rand*(gCaT_max-gCaT_min) + gCaT_min;
    g_CaS(ii) = rand*(gCaS_max-gCaS_min) + gCaS_min;
    g_A(ii) = rand*(gA_max-gA_min) + gA_min;
    g_KCa(ii) = rand*(gKCa_max-gKCa_min) + gKCa_min;
    g_K(ii) = rand*(gK_max-gK_min) + gK_min;
    g_H(ii) = rand*(gH_max-gH_min) + gH_min;
    
    figure;set(gcf,'position',[0 0 2000 400]);

    gr = 1;
    for gg=0:0.05:0.5
        [tot_T,V] = run_STG_model(g_Na(ii),g_K(ii),g_CaT(ii),g_CaS(ii),g_KCa(ii),g_A(ii),g_H(ii),gg,0); %identical neurons
        VV{ii,gr} = V;
        if (gr<11)
            subplot(2,5,gr);plot(tot_T,V(:,1),'k');
            hold on;plot(tot_T,V(:,2),'r');drawnow;hold off;
            xlim([1500 2000])
        end
        gr = gr+1;
    end
    %save('temp_spikers_onlygH.mat','VV','g_Na','g_K','g_CaT','g_CaS','g_A','g_H','g_KCa','-v7.3');
end

%temp_spikers.mat has [800,1200], [0,6], [0,12], [20,130], [100,500], [90,120], [1,1.5]