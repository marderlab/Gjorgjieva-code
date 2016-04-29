tt = 0:0.1:5000;
figure;
kk=1;

vec_test = [25256 35943 8498 33724 6904 13625];
vec_test_V = [879 1295];
vec_test2 = [14973 12392 24635 3575 19668 43457 59045 10685 40404 18745 51845 35478 21858 51564 15902 65410 23438 53461 8111 45887 ...
             64973 50771 5820 20572 5884 14921 49320 61157 24631 34513 47279 15786 21760 51391 42901 41354 41379 47513 ...
             8226 10132 10155 10198 12699 12710 12717 15970 27317 27375 30165 33093 36348 41908 42299]; 
ll=1; %1:6 8:9 11:14 16:21 24 26:29 33:34 36:45 47:53

for ii=pairs(record2(vec_test2(ll)),:)
    g_Na = round(g(ii,2)*10^6)/10^6;
    g_CaT = round(g(ii,3)*10^6)/10^6;
    g_CaS = round(g(ii,4)*10^6)/10^6;
    g_A = round(g(ii,5)*10^6)/10^6;
    g_KCa = round(g(ii,6)*10^6)/10^6;
    g_K = round(g(ii,7)*10^6)/10^6;
    g_H = round(g(ii,8)*10^6)/10^6;
    format short g
    gg = [g_Na g_CaT g_CaS g_A g_KCa g_K g_H];
    g_vec(kk,:) = gg;
    g_L = 0.01;
    fname = sprintf('V_gNa_%g_gCaT_%g_gCaS_%g_gA_%g_gKCa_%g_gK_%g_gH_%g_gL_%g.dat',g_Na,g_CaT,g_CaS,g_A,g_KCa,g_K,g_H,g_L);

    %cd data_files_from_cluster/
    cd data_files_highH
    try
        V = dlmread(fname);
        %plot(t,V); pause(0.1);
    catch
        try
            structur = dir(sprintf('V_gNa_%g_*.dat',g_Na));
            fname = structur.name;
            V = dlmread(fname);
            %plot(t,V); pause(0.1);
        catch
            try
                structur = dir(sprintf('V_*gCaS_%g_*.dat',g_CaS));
                fname = structur.name;
                V = dlmread(fname);
                %plot(t,V); pause(0.1);
            catch
                vec_weird = [vec_weird ii];
            end
        end
    end
    cd ..
    subplot(2,2,kk);plot(tt,V,'b');xlim([4500 5000]);ylim([-80 60])
    title(sprintf('pair %d (cell %d)',vec_test2(ll),pairs(record2(vec_test2(ll)),kk)))
    subplot(2,2,kk+2);bar([gg(1)/100 gg(2:3) gg(4)/10 gg(5)/20  gg(6)/10 gg(7)*5]);ylim([0 25]) %[0 14]
    set(gca,'xticklabel',{'g_{Na}/100','g_{CaT}','g_{CaS}','g_A/10','g_{KCa}/20','g_{K}/10','g_H*5'})
    kk = kk+1;
end
subplot(2,2,3);    
title(sprintf('MSE = %g',mean((g_vec(1,:)-g_vec(2,:)).^2)))
subplot(2,2,4);    
title(sprintf('L1 = %g',mean(abs(g_vec(1,:)-g_vec(2,:)))))
set(gcf,'position',[0 0 800 600]);
%%
%calculate random MSE
MSE = [];
L1 = [];
for kk=1:num_pairs
    ii=floor(rand*2000)+1;
    jj=floor(rand*2000)+1;
    g1 = g(ii,2:8); g2 = g(jj,2:8);
    MSE(kk) = mean((g1-g2).^2);
    L1(kk) = mean(abs(g1-g2));
end
