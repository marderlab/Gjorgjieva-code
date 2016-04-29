g = load('data_files_from_cluster/conductances_1_2000.dat');
%plot_grid(g(:,2:8),[800 1200; 0 6; 0 12; 20 130; 20 140; 90 120; 0.05 0.25]);

clear bp bp_temp
t=0:0.1:5000;

Num_sims = 2000;
%figure;
vec_weird = [];
weird = 0;
for ii=1:Num_sims
    if (mod(ii,100)==0)
        ii
    end
    %for some reason, %g did not work the same in C and Matlab so some
    %files won't load because last digit difference in the name
    g_Na = round(g(ii,2)*10^6)/10^6;
    g_CaT = round(g(ii,3)*10^6)/10^6;
    g_CaS = round(g(ii,4)*10^6)/10^6;
    g_A = round(g(ii,5)*10^6)/10^6;
    g_KCa = round(g(ii,6)*10^6)/10^6;
    g_K = round(g(ii,7)*10^6)/10^6;
    g_H = round(g(ii,8)*10^6)/10^6;
    g_L = 0.01;
    fname = sprintf('sp_gNa_%g_gCaT_%g_gCaS_%g_gA_%g_gKCa_%g_gK_%g_gH_%g_gL_%g.dat',g_Na,g_CaT,g_CaS,g_A,g_KCa,g_K,g_H,g_L);
    
    cd data_files_from_cluster/
    try
        sp = dlmread(fname);
    catch
        try
            structur = dir(sprintf('sp_gNa_%g_*.dat',g_Na));
            fname = structur.name;
            sp = dlmread(fname);
        catch
            try
                structur = dir(sprintf('sp_*gCaS_%g_*.dat',g_CaS));
                fname = structur.name;
                sp = dlmread(fname);
            catch
                vec_weird = [vec_weird ii];
                weird = 1;
            end
        end
    end
    cd ..
    if (weird == 0)
        bp_temp = burst_params(sp);
        bp(ii,:) = bp_temp{1};
    else
        bp(ii,:) = [NaN NaN NaN];
        weird = 0;
    end
end
%for some reason, before time 2000,000 seconds, this neuron is weird, but I
%don't have the voltage to prove it, look into it later. During the time
%for which voltage was recorded, the neuron is spiking at very high rate.
bp_temp = burst_params(sp(sp>2000000));
bp(632,:)=bp_temp;
%%

figure;subplot(3,1,1);hist(bp(:,1),50);tisize(bp)tle('# spikes/burst')
subplot(3,1,2);hist(bp(:,2),50);title('burst period (ms)')
subplot(3,1,3);hist(bp(:,3),50);title('duty cycle')
%%
vec_min_max = [];
for jj=1:size(bp,2)
    vec_min_max = [vec_min_max; min(bp(:,jj)) max(bp(:,jj))];
end

g_min_max = [800 1200; 0 6; 0 12; 20 130; 20 140; 90 120; 0.05 0.25];

A = g(:,2:8);
figure;
N = size(A,2);
for ii=1:N
    for jj=ii+1:N
        subplot(N-1,N-1,(ii-1)*N+jj - ii);scatter(A(:,ii),A(:,jj),10,bp(:,1),'filled');colorbar;
        caxis([min(bp(:,1)) quantile(bp(:,1),0.9)])
        xlim([g_min_max(ii,1) g_min_max(ii,2)])
        ylim([g_min_max(jj,1) g_min_max(jj,2)])
    end
end
subplot(6,6,1);xlabel('g_{Na}');
subplot(6,6,1);ylabel('g_{CaT}');
subplot(6,6,8);xlabel('g_{CaT}');
subplot(6,6,8);ylabel('g_{CaS}');
subplot(6,6,15);xlabel('g_{CaS}');
subplot(6,6,15);ylabel('g_{A}');
subplot(6,6,22);xlabel('g_{A}');
subplot(6,6,22);ylabel('g_{KCa}');
subplot(6,6,29);xlabel('g_{KCa}');
subplot(6,6,29);ylabel('g_{Kd}');
subplot(6,6,36);xlabel('g_{Kd}');
subplot(6,6,36);ylabel('g_{H}');
%%%%%%%%%%%%%%
figure;
N = size(A,2);
for ii=1:N
    for jj=ii+1:N
        subplot(N-1,N-1,(ii-1)*N+jj - ii);scatter(A(:,ii),A(:,jj),10,bp(:,2),'filled');colorbar
        caxis([min(bp(:,2)) quantile(bp(:,2),0.9)])
        xlim([g_min_max(ii,1) g_min_max(ii,2)])
        ylim([g_min_max(jj,1) g_min_max(jj,2)])
    end
end
subplot(6,6,1);xlabel('g_{Na}');
subplot(6,6,1);ylabel('g_{CaT}');
subplot(6,6,8);xlabel('g_{CaT}');
subplot(6,6,8);ylabel('g_{CaS}');
subplot(6,6,15);xlabel('g_{CaS}');
subplot(6,6,15);ylabel('g_{A}');
subplot(6,6,22);xlabel('g_{A}');
subplot(6,6,22);ylabel('g_{KCa}');
subplot(6,6,29);xlabel('g_{KCa}');
subplot(6,6,29);ylabel('g_{Kd}');
subplot(6,6,36);xlabel('g_{Kd}');
subplot(6,6,36);ylabel('g_{H}');
%%%%%%%%%%%%%%
figure;
N = size(A,2);
for ii=1:N
    for jj=ii+1:N
        subplot(N-1,N-1,(ii-1)*N+jj - ii);scatter(A(:,ii),A(:,jj),10,bp(:,3),'filled');colorbar
        caxis([min(bp(:,3)) quantile(bp(:,3),0.9)])
        xlim([g_min_max(ii,1) g_min_max(ii,2)])
        ylim([g_min_max(jj,1) g_min_max(jj,2)])
    end
end
subplot(6,6,1);xlabel('g_{Na}');
subplot(6,6,1);ylabel('g_{CaT}');
subplot(6,6,8);xlabel('g_{CaT}');
subplot(6,6,8);ylabel('g_{CaS}');
subplot(6,6,15);xlabel('g_{CaS}');
subplot(6,6,15);ylabel('g_{A}');
subplot(6,6,22);xlabel('g_{A}');
subplot(6,6,22);ylabel('g_{KCa}');
subplot(6,6,29);xlabel('g_{KCa}');
subplot(6,6,29);ylabel('g_{Kd}');
subplot(6,6,36);xlabel('g_{Kd}');
subplot(6,6,36);ylabel('g_{H}');
    