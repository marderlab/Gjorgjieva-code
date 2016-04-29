g = load('data_files_from_cluster/conductances_1_2000.dat');
%plot_grid(g(:,2:8),[800 1200; 0 6; 0 12; 20 130; 20 140; 90 120; 0.05 0.25]);

t=0:0.1:5000;

Num_sims = 2000;
%figure;
vec_weird = [];
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
    fname = sprintf('V_gNa_%g_gCaT_%g_gCaS_%g_gA_%g_gKCa_%g_gK_%g_gH_%g_gL_%g.dat',g_Na,g_CaT,g_CaS,g_A,g_KCa,g_K,g_H,g_L);
    
    cd data_files_from_cluster/
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
    most_hyper_V(ii) = min(V);    
end
figure;hist(most_hyper_V,50);

%%

A = g(:,2:8);
figure;
N = size(A,2);
for ii=1:N
    for jj=ii+1:N
        subplot(N-1,N-1,(ii-1)*N+jj - ii);scatter(A(:,ii),A(:,jj),10,most_hyper_V,'filled');colorbar;
        caxis([min(most_hyper_V) -58])
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
