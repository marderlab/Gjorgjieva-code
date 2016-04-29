function FR = FI_curves_corr_var(dt,sig)

%mu = [0.1:0.1:1.5];
mu = [0:0.05:1.5];
%sig = 0.0;
tc = 1.0;
  
gNa = 80;
rat1 = [0.05 0.1 0.2];
rat2 = [0.05 0.1 0.2 0.3];
rat3 = [0.00001 0.0001 0.001 0.0032 0.01];

FR = zeros(length(rat1), length(rat2), length(rat3), length(mu));
for jj=1:length(rat1)
    for kk=1:length(rat2)
        for tt=1:length(rat3)
            for ll=1:length(mu)
                [jj kk tt ll]
                %filename = sprintf('sig_0p5/spikes2_mean_%g_sig_%g_gNa_%g_rat1_%g_rat2_%g_rat3_%g_tc_%g_dt_%g.dat', mu(ll), sig, gNa, rat1(jj), rat2(kk), rat3(tt), tc, dt);
                filename = sprintf('/Volumes/Grass_Gjorgjieva/correl_transfer/sims_thr_detect/spikes_mean_%g_sig_%g_gNa_%g_rat1_%g_rat2_%g_rat3_%g_tc_%g_dt_%g.dat', mu(ll), sig, gNa, rat1(jj), rat2(kk), rat3(tt), tc, dt);                    
                try spikes = dlmread(filename);
                    FR(jj,kk,tt,ll) = length(spikes)/60;
                end
            end
        end
    end
end


