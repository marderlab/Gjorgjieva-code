tt = 0:0.1:5000;
figure;
kk=3;

vec_test = [25256 35943 8498 33724 6904 13625];
vec_test_V = [879 1295];
         
vec_test2 = [14973 12392 24635 3575 19668 43457 59045 10685 40404 18745 51845 35478 21858 51564 15902 65410 23438 53461 8111 45887 ...
             64973 50771 5820 20572 5884 14921 49320 61157 24631 34513 47279 15786 21760 51391 42901 41354 41379 47513 ... %spikers with different periods
             8226 10132 10155 10198 12699 12710 12717 15970 27317 27375 30165 33093 36348 41908 42299 ... %spikers with high g_H
             12570 12614 12622 17299 17316 17322 47830 51483 59185 59191 59196 59199 62986 62998 62999 63934 65091 65097 65866 ... %spikers with high g_A
             4993 5006 6386 6412 6426 45367 45377 45384 58368 58377 58384 59892 59899 62406 62419 63366 63375]; %spikers with low g_CaT and g_CaS and high g_KCa         
ll=12;
%ll=1: small regions where neuron 1 has bursts, but neuron 2 has very spread out bursts, almost like spikers - saved Voltage
%ll=2: only one point like ll=1, the rest is spikers

load pairs_Feb16.mat

format short g
tau = 100;
Esyn = -78;
Vhalf = 45;
gsyn1 = 0.4;
gsyn2 = 0.26;

vec = pairs(record2(vec_test2(ll)),:);
for ii=1:2
    g_Na(ii) = round(g(vec(ii),2)*10^6)/10^6;
    g_CaT(ii) = round(g(vec(ii),3)*10^6)/10^6;
    g_CaS(ii) = round(g(vec(ii),4)*10^6)/10^6;
    g_A(ii) = round(g(vec(ii),5)*10^6)/10^6;
    g_KCa(ii) = round(g(vec(ii),6)*10^6)/10^6;
    g_K(ii) = round(g(vec(ii),7)*10^6)/10^6;
    g_H(ii) = round(g(vec(ii),8)*10^6)/10^6;
end
%cd data_half_centers/
cd data_files_highH
fname = sprintf('V_gNa_%g_%g_gCaT_%g_%g_gCaS_%g_%g_gA_%g_%g_gKCa_%g_%g_gK_%g_%g_gH_%g_%g_gL_0.01_0.01_gsyn_%g_%g_Vhalf_%g_tsyn_%g_Esyn_%g.dat',g_Na(1),g_Na(2),g_CaT(1),g_CaT(2),g_CaS(1),g_CaS(2),g_A(1),g_A(2),g_KCa(1),g_KCa(2),g_K(1),g_K(2),g_H(1),g_H(2),gsyn1,gsyn2,Vhalf,tau,Esyn);


try
    V = dlmread(fname);
    %plot(t,V); pause(0.1);
catch
    try
        structur = dir(sprintf('V_gNa_%g_%g_*gsyn_%g_%g*.dat',g_Na(1),g_Na(2),gsyn1,gsyn2));
        fname = structur.name;
        V = dlmread(fname);
        %plot(t,V); pause(0.1);
    catch
        try
            structur = dir(sprintf('V_*gCaS_%g_%g_*gsyn_%g_%g*.dat',g_CaS(1),g_CaS(2),gsyn1,gsyn2));
            fname = structur.name;
            V = dlmread(fname);
            %plot(t,V); pause(0.1);
        catch
            try
                structur = dir(sprintf('V_*gA_%g_%g*gsyn_%g_%g*.dat',g_A(1),g_A(2),gsyn1,gsyn2));
                fname = structur.name;
                V = dlmread(fname);
                %plot(t,V); pause(0.1);
            catch
                vec_weird = [vec_weird ii];
            end
        end
    end
end
cd ..


V1 = V(V(:,1)==0,2);
V2 = V(V(:,1)==1,2);

subplot(2,1,1);
plot(tt,V1,'b','linewidth',2);xlim([4000 5000]);ylim([-80 60]);set(gca,'fontsize',20);
title(sprintf('pair %d (cell %d)',vec_test2(ll),pairs(record2(vec_test2(ll)),1)))
subplot(2,1,2);
hold on;plot(tt,V2,'r','linewidth',2);xlim([4000 5000]);ylim([-80 60])
title(sprintf('pair %d (cell %d)',vec_test2(ll),pairs(record2(vec_test2(ll)),2)))
%title(sprintf('pair %d (cells %d, %d)',vec_test2(ll),pairs(record2(vec_test2(ll)),1),pairs(record2(vec_test2(ll)),2)))
%set(gcf,'position',[0 0 1000 200]);
set(gca,'fontsize',20);set(gcf,'position',[0 0 1400 500]);