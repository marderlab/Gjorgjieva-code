%NOTE: BECAUSE THE CELLS ARE IDENTICAL/SYMMETRIC (INTRINSIC AND SYNAPTIC
%PROPERTIES) - THE PROPERTIES OF EACH CELL IN THE HALF-CENTER ARE IDENTICAL
%TO THE PROPERTIES OF THE OTHER CELL
sig = 1;
Tp = 5000; %10 sec, for the computation of the burst exclusion metric
thr_peak = 0.01;

t_vec = [1 10 100];
c_vec = [0:0.1:0.5 1];
%set of first two
for tt=1:length(t_vec)
    tau=t_vec(tt);
    for ccc=1:length(c_vec)
        corr = c_vec(ccc);
        cd ..
        name1 = sprintf('spikes1_mean_0_sig_%g_gNa_227.052_gCaT_2.7781_gCaS_3.3469_gA_30.4321_gKCa_121.124_gK_75.649_gH_1.3194_gL_0.1631_gsyn_0.005_0.005_corr_%g_tc_%d_dt_0.01_num_0.dat',sig,corr,tau);
        sp1 = dlmread(name1);
        name2 = sprintf('spikes2_mean_0_sig_%g_gNa_227.052_gCaT_2.7781_gCaS_3.3469_gA_30.4321_gKCa_121.124_gK_75.649_gH_1.3194_gL_0.1631_gsyn_0.005_0.005_corr_%g_tc_%d_dt_0.01_num_0.dat',sig,corr,tau);
        sp2 = dlmread(name2);
        cd data_base
        
        %[CV1, num_sp1, IBI1, burst_dur, first_spike, last_spike] = get_train_properties_sp(sp1,1);
        %[CV2, num_sp2, IBI2, burst_dur, first_spike, last_spike] = get_train_properties_sp(sp2,1);
        %figure(100);hold on;bar([ccc-0.2 ccc+0.2],[CV1 CV2],0.8)
        %figure(101);hold on;bar([ccc-0.2 ccc+0.2],[1/mean(IBI1) 1/mean(IBI2)],0.8)
        %figure(102);hold on;bar([ccc-0.2 ccc+0.2],[mean(num_sp1) mean(num_sp2)],0.8)
        %CV_pair1(tt,ccc) = CV1;
        %burst_freq_pair1(tt,ccc) = 1/mean(IBI1);
        %num_sp_pair1(tt,ccc) = mean(num_sp1);
        
        bin_width = 10; %ms, for binning spike trains to compute correlations
        offset = 0.5*10^4;
        time_vec = [0:bin_width:2000*10^3];
        T = 500; %ms (if bin_width = 10 ms, then T can be 500) OR (if bin_width = 7.5 ms, then T can be 750)
        corr_vec = [-T:bin_width:T];
        cell1_binned = histc(sp1-offset,time_vec); %
        cell2_binned = histc(sp2-offset,time_vec); %
        cc = xcov(cell1_binned,cell2_binned,T/bin_width,'coef');%/bin_width; %
        %figure;subplot(131);plot(corr_vec,cc/bin_width);xlim([-T T]);title(sprintf('tau=%g, corr=%g',tau,corr));
        
        %cc1 = xcov(cell1_binned,cell1_binned,T/bin_width,'coef'); %
        %subplot(132);plot(corr_vec,cc1,'r');xlim([-T T])
        %dc1 = diff(cc1);
        %note: to call something a peak, require that for two time points
        %before the derivative was positive, and for two time points after,
        %the derivative was negative
        %x1 = find(dc1(1:end-1)>thr_peak*abs(cc1(2:end-1)) | abs(dc1(2:end))>thr_peak*abs(cc1(2:end-1)))+1; 
        %x1 = intersect(x1,find(dc1(1:end-3)>0 & dc1(2:end-2)>0 & dc1(3:end-1)<0 & dc1(4:end)<0)+2);
        %coef1 = find(corr_vec(x1)>0);
        %temp = corr_vec(x1);
        %ibi_1 = temp(coef1(1))*0.001; %convert to seconds
        %figure;plot(corr_vec,cc1,'r');xlim([-T T]);hold on;plot(corr_vec(x1),cc1(x1),'o')

        %cc2 = xcov(cell2_binned,cell2_binned,T/bin_width,'coef'); %
        %subplot(133);plot(corr_vec,cc2,'k');xlim([-T T])
        %dc2 = diff(cc2);
        %x2 = find(dc2(1:end-1)>thr_peak*abs(cc2(2:end-1)) |  abs(dc2(2:end))>thr_peak*abs(cc2(2:end-1)))+1; 
        %x2 = intersect(x2,find(dc2(1:end-3)>0 & dc2(2:end-2)>0 & dc2(3:end-1)<0 & dc2(4:end)<0)+2);
        %coef2 = find(corr_vec(x2)>0);
        %temp = corr_vec(x2);
        %ibi_2 = temp(coef2(1))*0.001; %convert to seconds
        %figure;plot(corr_vec,cc2,'k');xlim([-T T]);hold on;plot(corr_vec(x2),cc2(x2),'o')
        
        %figure(103);hold on;bar([ccc-0.2 ccc+0.2],[1/ibi_1 1/ibi_2],0.8)
        %burst_freq2_pair1(tt,ccc) = 1/ibi_1;
        
        %%%now for the cross burst properties
        sign_at_0 = cc(T/bin_width+1); %value at 0
        %negative if anti-phase, otherwise positive
        %figure(104);hold on;bar(ccc,sign_at_0,0.8,'r')
        sign_at_0_pair1(tt,ccc) = sign_at_0;
        
        %dc = diff(cc);
        %x = find(dc(1:end-1)>thr_peak*abs(cc(2:end-1)) | abs(dc(2:end))>thr_peak*abs(cc(2:end-1)))+1; 
        %x = intersect(x,find(dc(1:end-3)>0 & dc(2:end-2)>0 & dc(3:end-1)<0 & dc(4:end)<0)+2);
        %figure;plot(corr_vec,cc,'b');xlim([-T T]);hold on;plot(corr_vec(x),cc(x),'o')
        %coef = find(corr_vec(x)>0);
        %temp = corr_vec(x);
        %temp_cc = cc(x);
        %te = 1;
        %to ensure that there weren't two inflections points near 0 - which
        %would be close in time, here consider the peak to be at least 
        %xx ms away from the center
        %the value of xx is kind of arbitrary - may need to also require
        %that the peak at this point is of a certain value (+ or -) but it
        %can be either + or - ...
        %while temp(coef(te))<100
        %    te = te+1;
        %end
        %peak1 = temp_cc(coef(te));
        %T1 = temp(coef(te));
        %
        %peak2 = temp_cc(coef(te+1));
        %T2 = temp(coef(te+1));
        %
        %figure(105);hold on;bar(ccc,peak1,0.8,'k')
        %figure(106);hold on;bar(ccc,T1,0.8,'c')
        
        %figure(107);hold on;bar(ccc,peak2,0.8,'k')
        %figure(108);hold on;bar(ccc,T2,0.8,'c')
        %peak1_pair1(tt,ccc) = peak1;
        %peak2_pair1(tt,ccc) = peak2;
        %time_to_burst_pair1(tt,ccc) = T1;
        
        chi_pair1(tt,ccc) = get_exclusion_metric(sp1,sp2,Tp);
        %prompt = 'Just enter a number to continue.';
        %x = input(prompt);
        %close all
    end
end

%%

%NOTE: BECAUSE THE CELLS ARE IDENTICAL/SYMMETRIC (INTRINSIC AND SYNAPTIC
%PROPERTIES) - THE PROPERTIES OF EACH CELL IN THE HALF-CENTER ARE IDENTICAL
%TO THE PROPERTIES OF THE OTHER CELL
%
%BUT ACTUALLY FOR THIS PARTICULAR SET THIS IS NOT THE CASE. I.E. SOMETIMES
%ONE OF THE CELLS BECOMES SILENT FOR LARGE PERIODS OF TIME
sig = 1;
Tp = 10000; %10 sec, for the computation of the burst exclusion metric
thr_peak = 0*0.01;

t_vec = [1 10 100];
c_vec = [0:0.1:0.5 1];

burst_freq_pair3 = zeros(length(t_vec),length(c_vec));
num_sp_pair3 = zeros(length(t_vec),length(c_vec));
burst_freq2_pair3 = zeros(length(t_vec),length(c_vec));
sign_at_0_pair3 = zeros(length(t_vec),length(c_vec));
peak1_pair3 = zeros(length(t_vec),length(c_vec));
time_to_burst_pair3 = zeros(length(t_vec),length(c_vec));
chi_pair3 = zeros(length(t_vec),length(c_vec));

%set of first two
for tt=1:length(t_vec)
    tau=t_vec(tt);
    for ccc=1:length(c_vec)
        corr = c_vec(ccc);
        run_further = 1;
        cd ..
        name1 = sprintf('spikes1_mean_0_sig_%g_gNa_227.052_gCaT_2.7781_gCaS_3.3469_gA_30.4321_gKCa_121.124_gK_75.649_gH_1.3194_gL_0.1631_gsyn_0.2_0.2_corr_%g_tc_%d_dt_0.01_num_0.dat',sig,corr,tau);
        try sp1 = dlmread(name1); catch run_further = 0; end
        name2 = sprintf('spikes2_mean_0_sig_%g_gNa_227.052_gCaT_2.7781_gCaS_3.3469_gA_30.4321_gKCa_121.124_gK_75.649_gH_1.3194_gL_0.1631_gsyn_0.2_0.2_corr_%g_tc_%d_dt_0.01_num_0.dat',sig,corr,tau);
        try sp2 = dlmread(name2); catch run_further = 0; end
        cd data_base
        
        if (run_further > 0)
        
            %[CV1, num_sp1, IBI1, burst_dur, first_spike, last_spike] = get_train_properties_sp(sp1,1);
            %[CV2, num_sp2, IBI2, burst_dur, first_spike, last_spike] = get_train_properties_sp(sp2,1);
            %figure(100);hold on;bar([ccc-0.2 ccc+0.2],[CV1 CV2],0.8)
            %figure(101);hold on;bar([ccc-0.2 ccc+0.2],[1/mean(IBI1) 1/mean(IBI2)],0.8)
            %figure(102);hold on;bar([ccc-0.2 ccc+0.2],[mean(num_sp1) mean(num_sp2)],0.8)
            %CV_pair3(tt,ccc) = CV1;
            %burst_freq_pair3(tt,ccc) = 1/mean(IBI1);
            %num_sp_pair3(tt,ccc) = mean(num_sp1);

            bin_width = 10; %ms, for binning spike trains to compute correlations
            offset = 0.5*10^4;
            time_vec = [0:bin_width:2000*10^3];
            T = 500; %ms (if bin_width = 10 ms, then T can be 500) OR (if bin_width = 7.5 ms, then T can be 750)
            corr_vec = [-T:bin_width:T];
            cell1_binned = histc(sp1-offset,time_vec); %
            cell2_binned = histc(sp2-offset,time_vec); %
            cc = xcov(cell1_binned,cell2_binned,T/bin_width,'coef');%/bin_width; %
            %figure;subplot(131);plot(corr_vec,cc/bin_width);xlim([-T T]);title(sprintf('tau=%g, corr=%g',tau,corr));

            %cc1 = xcov(cell1_binned,cell1_binned,T/bin_width,'coef'); %
            %subplot(132);plot(corr_vec,cc1,'r');xlim([-T T])
            %dc1 = diff(cc1);
            %note: to call somethign a peak, require that for two time points
            %before the derivative was positive, and for two time points after,
            %the derivative was negative
            %x1 = find(dc1(1:end-1)>thr_peak*abs(cc1(2:end-1)) | abs(dc1(2:end))>thr_peak*abs(cc1(2:end-1)))+1; 
            %x1 = intersect(x1,find(dc1(1:end-3)>0 & dc1(2:end-2)>0 & dc1(3:end-1)<0 & dc1(4:end)<0)+2);
            %coef1 = find(corr_vec(x1)>0);
            %temp = corr_vec(x1);
            %ibi_1 = temp(coef1(1))*0.001; %convert to seconds
            %figure;plot(corr_vec,cc1,'r');xlim([-T T]);hold on;plot(corr_vec(x1),cc1(x1),'o')

            %cc2 = xcov(cell2_binned,cell2_binned,T/bin_width,'coef'); %
            %subplot(133);plot(corr_vec,cc2,'k');xlim([-T T])
            %dc2 = diff(cc2);
            %x2 = find(dc2(1:end-1)>thr_peak*abs(cc2(2:end-1)) |  abs(dc2(2:end))>thr_peak*abs(cc2(2:end-1)))+1; 
            %x2 = intersect(x2,find(dc2(1:end-3)>0 & dc2(2:end-2)>0 & dc2(3:end-1)<0 & dc2(4:end)<0)+2);
            %coef2 = find(corr_vec(x2)>0);
            %temp = corr_vec(x2);
            %ibi_2 = temp(coef2(1))*0.001; %convert to seconds
            %figure;plot(corr_vec,cc2,'k');xlim([-T T]);hold on;plot(corr_vec(x2),cc2(x2),'o')

            %figure(103);hold on;bar([ccc-0.2 ccc+0.2],[1/ibi_1 1/ibi_2],0.8)
            %burst_freq2_pair3(tt,ccc) = 1/ibi_1;

            %%%now for the cross burst properties
            sign_at_0 = cc(T/bin_width+1); %value at 0
            %negative if anti-phase, otherwise positive
            %figure(104);hold on;bar(ccc,sign_at_0,0.8,'r')
            sign_at_0_pair3(tt,ccc) = sign_at_0;

            %dc = diff(cc);
            %x = find(dc(1:end-1)>thr_peak*abs(cc(2:end-1)) | abs(dc(2:end))>thr_peak*abs(cc(2:end-1)))+1; 
            %x = intersect(x,find(dc(1:end-3)>0 & dc(2:end-2)>0 & dc(3:end-1)<0 & dc(4:end)<0)+2);
            %figure;plot(corr_vec,cc,'b');xlim([-T T]);hold on;plot(corr_vec(x),cc(x),'o')
            %coef = find(corr_vec(x)>0);
            %temp = corr_vec(x);
            %temp_cc = cc(x);
            %te = 1;
            %to ensure that there weren't two inflections points near 0 - which
            %would be close in time, here consider the peak to be at least 
            %xx ms away from the center
            %the value of xx is kind of arbitrary - may need to also require
            %that the peak at this point is of a certain value (+ or -) but it
            %can be either + or - ...
            %while temp(coef(te))<100
            %    te = te+1;
            %end
            %peak1 = temp_cc(coef(te));
            %T1 = temp(coef(te));
            %
            %peak2 = temp_cc(coef(te+1));
            %T2 = temp(coef(te+1));
            %
            %figure(105);hold on;bar(ccc,peak1,0.8,'k')
            %figure(106);hold on;bar(ccc,T1,0.8,'c')

            %figure(107);hold on;bar(ccc,peak2,0.8,'k')
            %figure(108);hold on;bar(ccc,T2,0.8,'c')
            %peak1_pair3(tt,ccc) = peak1;
            %peak2_pair3(tt,ccc) = peak2;
            %time_to_burst_pair3(tt,ccc) = T1;

            chi_pair3(tt,ccc) = get_exclusion_metric(sp1,sp2,Tp);
            %prompt = 'Just enter a number to continue.';
            %x = input(prompt);
            %close all
        end
    end
end
%%
sig = 1;
Tp = 5000; %10 sec, for the computation of the burst exclusion metric
thr_peak = 0.01;

t_vec = [1 10 100];
c_vec = [0:0.1:0.5 1];
thr_peak = 0.01;
for tt=1:length(t_vec)
    tau=t_vec(tt);
    for ccc=1:length(c_vec)
        corr = c_vec(ccc);
        cd ..
        name1 = sprintf('spikes1_mean_5_sig_%g_gNa_227.052_gCaT_2.7781_gCaS_3.3469_gA_45_gKCa_121.124_gK_75.649_gH_1.3194_gL_0.1631_gsyn_0.2_0.2_corr_%g_tc_%d_dt_0.01_num_0.dat',sig,corr,tau);
        sp1 = dlmread(name1);
        name2 = sprintf('spikes2_mean_5_sig_%g_gNa_227.052_gCaT_2.7781_gCaS_3.3469_gA_45_gKCa_121.124_gK_75.649_gH_1.3194_gL_0.1631_gsyn_0.2_0.2_corr_%g_tc_%d_dt_0.01_num_0.dat',sig,corr,tau);
        sp2 = dlmread(name2);
        cd data_base
        
        %[CV1, num_sp1, IBI1, burst_dur, first_spike, last_spike] = get_train_properties_sp(sp1,1);
        %[CV2, num_sp2, IBI2, burst_dur, first_spike, last_spike] = get_train_properties_sp(sp2,1);
        %figure(200);hold on;bar([ccc-0.2 ccc+0.2],[CV1 CV2],0.8)
        %figure(201);hold on;bar([ccc-0.2 ccc+0.2],[1/mean(IBI1) 1/mean(IBI2)],0.8)
        %figure(202);hold on;bar([ccc-0.2 ccc+0.2],[mean(num_sp1) mean(num_sp2)],0.8)
        %CV_pair2(tt,ccc) = CV1;
        %burst_freq_pair2(tt,ccc) = 1/mean(IBI1);
        %num_sp_pair2(tt,ccc) = mean(num_sp1);
        
        bin_width = 10; %ms, for binning spike trains to compute correlations
        offset = 0.5*10^4;
        time_vec = [0:bin_width:2000*10^3];
        T = 500; %ms (if bin_width = 10 ms, then T can be 500) OR (if bin_width = 7.5 ms, then T can be 750)
        corr_vec = [-T:bin_width:T];
        cell1_binned = histc(sp1-offset,time_vec); %
        cell2_binned = histc(sp2-offset,time_vec); %
        cc = xcov(cell1_binned,cell2_binned,T/bin_width,'coef');%/bin_width; %
        %figure;subplot(131);plot(corr_vec,cc/bin_width);xlim([-T T]);title(sprintf('tau=%g, corr=%g',tau,corr));
        
        %cc1 = xcov(cell1_binned,cell1_binned,T/bin_width,'coef'); %
        %subplot(132);plot(corr_vec,cc1,'r');xlim([-T T])
        %dc1 = diff(cc1);
        %note: to call somethign a peak, require that for two time points
        %before the derivative was positive, and for two time points after,
        %the derivative was negative
        %x1 = find(dc1(1:end-1)>thr_peak*abs(cc1(2:end-1)) | abs(dc1(2:end))>thr_peak*abs(cc1(2:end-1)))+1; 
        %x1 = intersect(x1,find(dc1(1:end-3)>0 & dc1(2:end-2)>0 & dc1(3:end-1)<0 & dc1(4:end)<0)+2);
        %coef1 = find(corr_vec(x1)>0);
        %temp = corr_vec(x1);
        %ibi_1 = temp(coef1(1))*0.001; %convert to seconds
        %figure;plot(corr_vec,cc1,'r');xlim([-T T]);hold on;plot(corr_vec(x1),cc1(x1),'o')

        %cc2 = xcov(cell2_binned,cell2_binned,T/bin_width,'coef'); %
        %subplot(133);plot(corr_vec,cc2,'k');xlim([-T T])
        %dc2 = diff(cc2);
        %x2 = find(dc2(1:end-1)>thr_peak*abs(cc2(2:end-1)) |  abs(dc2(2:end))>thr_peak*abs(cc2(2:end-1)))+1; 
        %x2 = intersect(x2,find(dc2(1:end-3)>0 & dc2(2:end-2)>0 & dc2(3:end-1)<0 & dc2(4:end)<0)+2);
        %coef2 = find(corr_vec(x2)>0);
        %temp = corr_vec(x2);
        %ibi_2 = temp(coef2(1))*0.001; %convert to seconds
        %figure;plot(corr_vec,cc2,'k');xlim([-T T]);hold on;plot(corr_vec(x2),cc2(x2),'o')
        
        %figure(203);hold on;bar([ccc-0.2 ccc+0.2],[1/ibi_1 1/ibi_2],0.8)
        %burst_freq2_pair2(tt,ccc) = 1/ibi_1;
        
        %%%now for the cross burst properties
        sign_at_0 = cc(T/bin_width+1); %value at 0
        %negative if anti-phase, otherwise positive
        %figure(204);hold on;bar(ccc,sign_at_0,0.8,'r')
        sign_at_0_pair2(tt,ccc) = sign_at_0;
        
        %dc = diff(cc);
        %x = find(dc(1:end-1)>thr_peak*abs(cc(2:end-1)) | abs(dc(2:end))>thr_peak*abs(cc(2:end-1)))+1; 
        %x = intersect(x,find(dc(1:end-3)>0 & dc(2:end-2)>0 & dc(3:end-1)<0 & dc(4:end)<0)+2);
        %figure;plot(corr_vec,cc,'b');xlim([-T T]);hold on;plot(corr_vec(x),cc(x),'o')
        %coef = find(corr_vec(x)>0);
        %temp = corr_vec(x);
        %temp_cc = cc(x);
        %te = 1;
        %to ensure that there weren't two inflections points near 0 - which
        %would be close in time, here consider the peak to be at least 
        %xx ms away from the center
        %the value of xx is kind of arbitrary - may need to also require
        %that the peak at this point is of a certain value (+ or -) but it
        %can be either + or - ...
        %while temp(coef(te))<100
        %    te = te+1;
        %end
        %peak1 = temp_cc(coef(te));
        %T1 = temp(coef(te));
        %
        %peak2 = temp_cc(coef(te+1));
        %T2 = temp(coef(te+1));
        %
        %figure(205);hold on;bar(ccc,peak1,0.8,'k')
        %figure(206);hold on;bar(ccc,T1,0.8,'c')
        
        %figure(207);hold on;bar(ccc,peak2,0.8,'k')
        %figure(208);hold on;bar(ccc,T2,0.8,'c')
        %peak1_pair2(tt,ccc) = peak1;
        %peak2_pair2(tt,ccc) = peak2;
        %time_to_burst_pair2(tt,ccc) = T1;
        
        chi_pair2(tt,ccc) = get_exclusion_metric(sp1,sp2,Tp);
        %prompt = 'Just enter a number to continue.';
        %x = input(prompt);
        %close all
    end
end


% %%
% for tau=[1 10 100]
%     for corr=[0 0.5 1]
%         name1 = sprintf('spikes1_mean_5_sig_%g_gNa_227.052_gCaT_2.7781_gCaS_3.3469_gA_45_gKCa_121.124_gK_75.649_gH_1.3194_gL_0.1631_gsyn_0.2_0.2_corr_%g_tc_%d_dt_0.01_num_0.dat',sig,corr,tau);
%         sp1 = dlmread(name1);
%         name2 = sprintf('spikes2_mean_5_sig_%g_gNa_227.052_gCaT_2.7781_gCaS_3.3469_gA_45_gKCa_121.124_gK_75.649_gH_1.3194_gL_0.1631_gsyn_0.2_0.2_corr_%g_tc_%d_dt_0.01_num_0.dat',sig,corr,tau);
%         sp2 = dlmread(name2);
% 
%         bin_width = 5; %ms, for binning spike trains to compute correlations
%         offset = 0.5*10^4;
%         time_vec = [0:bin_width:2000*10^3];
%         T = 500; %s
%         corr_vec = [-T:bin_width:T];
%         cell1_binned = histc(sp1-offset,time_vec); %
%         cell2_binned = histc(sp2-offset,time_vec); %
%         cc = xcov(cell1_binned,cell2_binned,T/bin_width,'coef'); %
%         figure;subplot(131);plot(corr_vec,cc/bin_width);xlim([-T T]);title(sprintf('tau=%g, corr=%g',tau,corr));       
% 
%         cc1 = xcov(cell1_binned,cell1_binned,T/bin_width,'coef'); %
%         subplot(132);plot(corr_vec,cc1,'r');xlim([-T T])
% 
%         cc2 = xcov(cell2_binned,cell2_binned,T/bin_width,'coef'); %
%         subplot(133);plot(corr_vec,cc2,'k');xlim([-T T])
%     end
% end