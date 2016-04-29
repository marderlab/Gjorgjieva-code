sig = 1;

%set of first two
for tau=[1 10 100]
    for corr=[0 0.5 1]
        name1 = sprintf('spikes1_mean_0_sig_%g_gNa_227.052_gCaT_2.7781_gCaS_3.3469_gA_30.4321_gKCa_121.124_gK_75.649_gH_1.3194_gL_0.1631_gsyn_0.005_0.005_corr_%g_tc_%d_dt_0.01_num_0.dat',sig,corr,tau);
        sp1 = dlmread(name1);
        name2 = sprintf('spikes2_mean_0_sig_%g_gNa_227.052_gCaT_2.7781_gCaS_3.3469_gA_30.4321_gKCa_121.124_gK_75.649_gH_1.3194_gL_0.1631_gsyn_0.005_0.005_corr_%g_tc_%d_dt_0.01_num_0.dat',sig,corr,tau);
        sp2 = dlmread(name2);

        bin_width = 1; %ms, for binning spike trains to compute correlations
        offset = 0.5*10^4;
        time_vec = [0:bin_width:2000*10^3];
        T = 500; %s
        corr_vec = [-T:bin_width:T];
        cell1_binned = histc(sp1-offset,time_vec); %
        cell2_binned = histc(sp2-offset,time_vec); %
        cc = xcov(cell1_binned,cell2_binned,T/bin_width,'coef'); %
        figure;subplot(131);plot(corr_vec,cc);xlim([-T T]);title(sprintf('tau=%g, corr=%g',tau,corr));
        %
        bin_width = 10; %ms, for binning spike trains to compute correlations
        offset = 0.5*10^4;
        time_vec = [0:bin_width:2000*10^3];
        T = 500; %s
        corr_vec = [-T:bin_width:T];
        cell1_binned = histc(sp1-offset,time_vec); %
        cell2_binned = histc(sp2-offset,time_vec); %
        cc = xcov(cell1_binned,cell2_binned,T/bin_width,'coef'); %
        hold on;plot(corr_vec,cc/10,'c');xlim([-T T])

        cc1 = xcov(cell1_binned,cell1_binned,T/bin_width,'coef'); %
        subplot(132);plot(corr_vec,cc1,'r');xlim([-T T])

        cc2 = xcov(cell2_binned,cell2_binned,T/bin_width,'coef'); %
        subplot(133);plot(corr_vec,cc2,'k');xlim([-T T])
    end
end
%%
for tau=[1 10 100]
    for corr=[0 0.5 1]
        name1 = sprintf('spikes1_mean_5_sig_%g_gNa_227.052_gCaT_2.7781_gCaS_3.3469_gA_45_gKCa_121.124_gK_75.649_gH_1.3194_gL_0.1631_gsyn_0.2_0.2_corr_%g_tc_%d_dt_0.01_num_0.dat',sig,corr,tau);
        sp1 = dlmread(name1);
        name2 = sprintf('spikes2_mean_5_sig_%g_gNa_227.052_gCaT_2.7781_gCaS_3.3469_gA_45_gKCa_121.124_gK_75.649_gH_1.3194_gL_0.1631_gsyn_0.2_0.2_corr_%g_tc_%d_dt_0.01_num_0.dat',sig,corr,tau);
        sp2 = dlmread(name2);

        bin_width = 1; %ms, for binning spike trains to compute correlations
        offset = 0.5*10^4;
        time_vec = [0:bin_width:2000*10^3];
        T = 500; %s
        corr_vec = [-T:bin_width:T];
        cell1_binned = histc(sp1-offset,time_vec); %
        cell2_binned = histc(sp2-offset,time_vec); %
        cc = xcov(cell1_binned,cell2_binned,T/bin_width,'coef'); %
        figure;subplot(131);plot(corr_vec,cc);xlim([-T T]);title(sprintf('tau=%g, corr=%g',tau,corr));
        %
        bin_width = 10; %ms, for binning spike trains to compute correlations
        offset = 0.5*10^4;
        time_vec = [0:bin_width:2000*10^3];
        T = 500; %s
        corr_vec = [-T:bin_width:T];
        cell1_binned = histc(sp1-offset,time_vec); %
        cell2_binned = histc(sp2-offset,time_vec); %
        cc = xcov(cell1_binned,cell2_binned,T/bin_width,'coef'); %
        hold on;plot(corr_vec,cc/10,'r');xlim([-T T])

        cc1 = xcov(cell1_binned,cell1_binned,T/bin_width,'coef'); %
        subplot(132);plot(corr_vec,cc1,'r');xlim([-T T])

        cc2 = xcov(cell2_binned,cell2_binned,T/bin_width,'coef'); %
        subplot(133);plot(corr_vec,cc2,'k');xlim([-T T])
    end
end
%%
%set of second two
for tau=[1 10 100]
    for corr=[0 0.5 1]
        name1 = sprintf('spikes1_mean_1.1_sig_%g_gNa_227.052_gCaT_2.7781_gCaS_3.3469_gA_30.4321_gKCa_121.124_gK_75.649_gH_1.3194_gL_0.1631_gsyn_0.01_0.01_corr_%g_tc_%d_dt_0.01_num_0.dat',sig,corr,tau);
        sp1 = dlmread(name1);
        name2 = sprintf('spikes2_mean_1.1_sig_%g_gNa_227.052_gCaT_2.7781_gCaS_3.3469_gA_30.4321_gKCa_121.124_gK_75.649_gH_1.3194_gL_0.1631_gsyn_0.01_0.01_corr_%g_tc_%d_dt_0.01_num_0.dat',sig,corr,tau);
        sp2 = dlmread(name2);

        bin_width = 1; %ms, for binning spike trains to compute correlations
        offset = 0.5*10^4;
        time_vec = [0:bin_width:2000*10^3];
        T = 500; %s
        corr_vec = [-T:bin_width:T];
        cell1_binned = histc(sp1-offset,time_vec); %
        cell2_binned = histc(sp2-offset,time_vec); %
        cc = xcov(cell1_binned,cell2_binned,T/bin_width,'coef'); %
        figure;subplot(131);plot(corr_vec,cc);xlim([-T T]);title(sprintf('tau=%g, corr=%g',tau,corr));
        %
        bin_width = 10; %ms, for binning spike trains to compute correlations
        offset = 0.5*10^4;
        time_vec = [0:bin_width:2000*10^3];
        T = 500; %s
        corr_vec = [-T:bin_width:T];
        cell1_binned = histc(sp1-offset,time_vec); %
        cell2_binned = histc(sp2-offset,time_vec); %
        cc = xcov(cell1_binned,cell2_binned,T/bin_width,'coef'); %
        hold on;plot(corr_vec,cc/10,'c');xlim([-T T])

        cc1 = xcov(cell1_binned,cell1_binned,T/bin_width,'coef'); %
        subplot(132);plot(corr_vec,cc1,'r');xlim([-T T])

        cc2 = xcov(cell2_binned,cell2_binned,T/bin_width,'coef'); %
        subplot(133);plot(corr_vec,cc2,'k');xlim([-T T])
    end
end

%%
for tau=[1 10 100]
    for corr=[0 0.5 1]
        name1 = sprintf('spikes1_mean_0_sig_%g_gNa_227.052_gCaT_2.2781_gCaS_2.8469_gA_30.4321_gKCa_1121.12_gK_75.649_gH_1.3194_gL_0.1631_gsyn_0.01_0.01_corr_%g_tc_%d_dt_0.01_num_0.dat',sig,corr,tau);
        sp1 = dlmread(name1);
        name2 = sprintf('spikes2_mean_0_sig_%g_gNa_227.052_gCaT_2.2781_gCaS_2.8469_gA_30.4321_gKCa_1121.12_gK_75.649_gH_1.3194_gL_0.1631_gsyn_0.01_0.01_corr_%g_tc_%d_dt_0.01_num_0.dat',sig,corr,tau);
        sp2 = dlmread(name2);

        bin_width = 1; %ms, for binning spike trains to compute correlations
        offset = 0.5*10^4;
        time_vec = [0:bin_width:2000*10^3];
        T = 500; %s
        corr_vec = [-T:bin_width:T];
        cell1_binned = histc(sp1-offset,time_vec); %
        cell2_binned = histc(sp2-offset,time_vec); %
        cc = xcov(cell1_binned,cell2_binned,T/bin_width,'coef'); %
        figure;subplot(131);plot(corr_vec,cc);xlim([-T T]);title(sprintf('tau=%g, corr=%g',tau,corr));
        %
        bin_width = 10; %ms, for binning spike trains to compute correlations
        offset = 0.5*10^4;
        time_vec = [0:bin_width:2000*10^3];
        T = 500; %s
        corr_vec = [-T:bin_width:T];
        cell1_binned = histc(sp1-offset,time_vec); %
        cell2_binned = histc(sp2-offset,time_vec); %
        cc = xcov(cell1_binned,cell2_binned,T/bin_width,'coef'); %
        hold on;plot(corr_vec,cc/10,'c');xlim([-T T])

        cc1 = xcov(cell1_binned,cell1_binned,T/bin_width,'coef'); %
        subplot(132);plot(corr_vec,cc1,'r');xlim([-T T])

        cc2 = xcov(cell2_binned,cell2_binned,T/bin_width,'coef'); %
        subplot(133);plot(corr_vec,cc2,'k');xlim([-T T])
    end
end

%%
%set of second two
for tau=[1 10 100]
    for corr=[0 0.5 1]
        name1 = sprintf('spikes1_mean_3_sig_%g_gNa_227.052_gCaT_2.7781_gCaS_3.3469_gA_30.4321_gKCa_121.124_gK_75.649_gH_1.3194_gL_0.1631_gsyn_0.01_0.01_corr_%g_tc_%d_dt_0.01_num_0.dat',sig,corr,tau);
        sp1 = dlmread(name1);
        name2 = sprintf('spikes2_mean_3_sig_%g_gNa_227.052_gCaT_2.7781_gCaS_3.3469_gA_30.4321_gKCa_121.124_gK_75.649_gH_1.3194_gL_0.1631_gsyn_0.01_0.01_corr_%g_tc_%d_dt_0.01_num_0.dat',sig,corr,tau);
        sp2 = dlmread(name2);

        bin_width = 1; %ms, for binning spike trains to compute correlations
        offset = 0.5*10^4;
        time_vec = [0:bin_width:2000*10^3];
        T = 500; %s
        corr_vec = [-T:bin_width:T];
        cell1_binned = histc(sp1-offset,time_vec); %
        cell2_binned = histc(sp2-offset,time_vec); %
        cc = xcov(cell1_binned,cell2_binned,T/bin_width,'coef'); %
        figure;subplot(131);plot(corr_vec,cc);xlim([-T T]);title(sprintf('tau=%g, corr=%g',tau,corr));
        %
        bin_width = 10; %ms, for binning spike trains to compute correlations
        offset = 0.5*10^4;
        time_vec = [0:bin_width:2000*10^3];
        T = 500; %s
        corr_vec = [-T:bin_width:T];
        cell1_binned = histc(sp1-offset,time_vec); %
        cell2_binned = histc(sp2-offset,time_vec); %
        cc = xcov(cell1_binned,cell2_binned,T/bin_width,'coef'); %
        hold on;plot(corr_vec,cc/10,'c');xlim([-T T])

        cc1 = xcov(cell1_binned,cell1_binned,T/bin_width,'coef'); %
        subplot(132);plot(corr_vec,cc1,'r');xlim([-T T])

        cc2 = xcov(cell2_binned,cell2_binned,T/bin_width,'coef'); %
        subplot(133);plot(corr_vec,cc2,'k');xlim([-T T])
    end
end
