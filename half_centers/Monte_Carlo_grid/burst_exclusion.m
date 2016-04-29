function chi = burst_exclusion(sp1,sp2,T_trial)
%function chi = burst_exclusion(V1,V2,t1)
%
%V1=dlmread('V1_mean_0_sig_0_gNa_227.052_gCaT_2.7781_gCaS_3.3469_gA_10_gKCa_121.124_gK_75.649_gH_1.3194_gL_0.1631_gsyn_0.001_0.001_corr_0_tc_0_dt_0.01_num_0.dat');
%V2=dlmread('V2_mean_0_sig_0_gNa_227.052_gCaT_2.7781_gCaS_3.3469_gA_10_gKCa_121.124_gK_75.649_gH_1.3194_gL_0.1631_gsyn_0.001_0.001_corr_0_tc_0_dt_0.01_num_0.dat');
%[~, ~, ~, ~, first_spike1, last_spike1] = get_train_properties(V1,t1);
%[~, ~, ~, ~, first_spike2, last_spike2] = get_train_properties(V2,t1);
[~, ~, ~, ~, first_spike1, last_spike1] = get_train_properties_sp(sp1,0); %the last parameter 0 is to not print the histogram of CVs
[~, ~, ~, ~, first_spike2, last_spike2] = get_train_properties_sp(sp2,0);

%T_trial = t1(end);
%T_trial = 5000;
%T_trial = 65000; %total duration of trial
Onet = 0; %the total time of activation of both cells
for ii=1:length(first_spike1)
    s1 = first_spike1(ii);
    e1 = last_spike1(ii);
    for jj=1:length(first_spike2)
        s2 = first_spike2(jj);
        e2 = last_spike2(jj);
        overlap = min(e1,e2)-max(s1,s2);
        Onet = Onet + overlap.*(overlap>0);
    end
end

t_1 = sum(last_spike1 - first_spike1); %total time of activity of cell 1
t_2 = sum(last_spike2 - first_spike2); %total time of activity of cell 2

if (t_1+t_2>T_trial)
    %overlap time if random
    Orand = min(t_1,t_2) - 0.5*(T_trial-max(t_1,t_2));
    %minimum possible overlap time
    Omin = (T_trial - t_1 - t_2).*(t_1+t_2>T_trial);
else
    Orand = min(t_1,t_2)^2/(2*(T_trial-max(t_1,t_2)));
    Omin = 0;
end
Omax = min(t_1,t_2); %max possible overlap time

if (Onet>Orand)
    chi = (Orand-Onet)/(Omax-Orand);
else
    chi = (Orand-Onet)/(Orand-Omin);
end



        