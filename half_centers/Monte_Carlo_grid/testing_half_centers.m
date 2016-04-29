[tot_T,V] = run_STG_model(800,90,2.2781,2.8469,800,30.4321,1.3194,0.4);set(gcf,'position',[0 0 1000 200]); %two sp bursts, spikers when uncoupled
[tot_T,V] = run_STG_model(800,90,2.2781,2.8469,700,30.4321,1.3194,0.4);set(gcf,'position',[0 0 1000 200]);
[tot_T,V] = run_STG_model(800,90,2.2781,2.8469,600,30.4321,1.3194,0.4);set(gcf,'position',[0 0 1000 200]);
[tot_T,V] = run_STG_model(800,90,2.2781,2.8469,500,30.4321,1.3194,0.35);set(gcf,'position',[0 0 1000 200]); %but not with 0.4 syn strength
[tot_T,V] = run_STG_model(800,90,2.2781,2.8469,400,30.4321,1.3194,0.35);set(gcf,'position',[0 0 1000 200]);
[tot_T,V] = run_STG_model(800,90,2.2781,2.8469,300,30.4321,1.3194,0.35);set(gcf,'position',[0 0 1000 200]); %still 2 sp bursts, single spikers
[tot_T,V] = run_STG_model(800,90,2.2781,2.8469,200,30.4321,1.3194,0.3);set(gcf,'position',[0 0 1000 200]); %still 2 sp bursts, single spikers

[tot_T,V] = run_STG_model(800,90,2.2781,2.8469,140,30.4321,1.3194,0.3);set(gcf,'position',[0 0 1000 200]); %even an occasional 3 spiker
[tot_T,V] = run_STG_model(800,90,2.2781,2.8469,140,30.4321,1.3194,0.31);set(gcf,'position',[0 0 1000 200]); %ones goes silent
%maybe it's hard to get bursters from spikers that are very different?
[tot_T,V] = run_STG_model(800,90,2.2781,2.8469,140,30.4321,1.3194,0.25);set(gcf,'position',[0 0 1000 200]); %2 sp, but also single neurons are 2 sp


%%
figure;
gr = 1;
for gg=0.:0.05:0.45 %g_Na,g_K,g_CaT,g_CaS,g_KCa,g_A,g_H
    [tot_T,V] = run_STG_model(800,100,2.5,6,130,90,1,gg,0); %identical neurons
    if (gr<11)
        subplot(2,5,gr);plot(tot_T,V(:,1),'k');
        hold on;plot(tot_T,V(:,2),'r');xlim([1500 2000]);drawnow;hold off;        
    end
    gr = gr+1;
end