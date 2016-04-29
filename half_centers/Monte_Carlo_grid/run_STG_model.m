function [tot_T,V] = run_STG_model(g_Na,g_K,g_CaT,g_CaS,g_KCa,g_A,g_H,gsyn,gapp,Vsyn,plot_traces)
%argh in the original runs - the excitatory input also meant excitatory
%synapses between the two neurons!
%%% Parameters
C = 1; 
VNa = 50; VK = -80; Vh = -20; Vleak = -50; %Vsyn = -78;%mV
gleak = 0.01; %mS/cm2
gNa = g_Na; gCaT = g_CaT; gCaS = g_CaS; 
gKd = g_K; gKCa = g_KCa; gA = g_A;  gH = g_H;
tausyn = 100;
R_F = 8.6174e-005; %needed for the computation of the revsersal for Ca
T = 10; tauCa = 20;

%Iapp = 0;

%%% Initial values
V(1,1:2) = [-50 -50];% + rand(1,2);
Ca(1,1:2) = 0.2*[1 1];
VCa = 120*ones(1,2); 
ICa = zeros(1,2);
mNa = Xinf(V(1,:),25.5,-5.29);
hNa = Xinf(V(1,:),48.9,5.18);
mCaT = Xinf(V(1,:),27.1,-7.2);
hCaT = Xinf(V(1,:),32.1,5.5);
mCaS = Xinf(V(1,:),33,-8.1);
hCaS = Xinf(V(1,:),60,6.2);
mA = Xinf(V(1,:),27.2,-8.7);
hA = Xinf(V(1,:),56.9,4.9);
mH = Xinf(V(1,:),70,6);
mKCa = m_KCainf(V(1,:),Ca(1,:));
mKd = Xinf(V(1,:),12.3,-11.8);
msyn = Xinf(V(1,:),45,-2);
Vsyn_HCO = -78;

%%% Simulation
tinit = 0;
tfinal = 10000; %ms
dt = 0.00125;
tot_T = tinit:dt:tfinal;

for tt=1:length(tot_T)
    VCa = 500.0*R_F*(T + 273.15)*log(3000./Ca(tt,:));
    
    Ca_inf = 0.05 - 0.94*ICa;
    Ca(tt+1,:) = Ca(tt,:) + (1-exp(-dt./tauCa)).*(Ca_inf - Ca(tt,:));
     
    mNainf = Xinf(V(tt,:),25.5,-5.29); taumNa = tauX(V(tt,:),1.32,1.26,120,-25);
    hNainf = Xinf(V(tt,:),48.9,5.18); tauhNa = tau_hNa(V(tt,:));
    mCaTinf = Xinf(V(tt,:),27.1,-7.2); taumCaT = tauX(V(tt,:),21.7,21.3,68.1,-20.5);
    hCaTinf = Xinf(V(tt,:),32.1,5.5); tauhCaT = tauX(V(tt,:),105,89.8,55,-16.9);
    mCaSinf = Xinf(V(tt,:),33,-8.1); taumCaS = tau_mCaS(V(tt,:));
    hCaSinf = Xinf(V(tt,:),60,6.2); tauhCaS = tau_hCaS(V(tt,:));
    mAinf = Xinf(V(tt,:),27.2,-8.7); taumA = tauX(V(tt,:),11.6,10.4,32.9,-15.2);
    hAinf = Xinf(V(tt,:),56.9,4.9); tauhA = tauX(V(tt,:),38.6,29.2,38.9,-26.5);
    mKCainf = m_KCainf(V(tt,:),Ca(tt,:)); taumKCa = tauX(V(tt,:),90.3,75.1,46,-22.7);
    mKdinf = Xinf(V(tt,:),12.3,-11.8); taumKd = tauX(V(tt,:),7.2,6.4,28.3,-19.2);
    mHinf = Xinf(V(tt,:),70,6); taumH = tauX(V(tt,:),272,-1499,42.2,-8.73);
    msyninf = Xinf(V(tt,:),45,-2);
    
    mNa = mNa + (1-exp(-dt./taumNa)).*(mNainf - mNa);
    hNa = hNa + (1-exp(-dt./tauhNa)).*(hNainf - hNa);
    mCaT = mCaT + (1-exp(-dt./taumCaT)).*(mCaTinf - mCaT);
    hCaT = hCaT + (1-exp(-dt./tauhCaT)).*(hCaTinf - hCaT);
    mCaS = mCaS + (1-exp(-dt./taumCaS)).*(mCaSinf - mCaS);
    hCaS = hCaS + (1-exp(-dt./tauhCaS)).*(hCaSinf - hCaS);
    mA = mA + (1-exp(-dt./taumA)).*(mAinf - mA);
    hA = hA + (1-exp(-dt./tauhA)).*(hAinf - hA); 
    mH = mH + (1-exp(-dt./taumH)).*(mHinf - mH);
    mKCa = mKCa + (1-exp(-dt./taumKCa)).*(mKCainf - mKCa);
    mKd = mKd + (1-exp(-dt./taumKd)).*(mKdinf - mKd);
    msyn = msyn + (1-exp(-dt./tausyn)).*(msyninf - msyn);
        
    sum_g = -(-gapp(tt,:)-gNa.*mNa.^3.*hNa -gCaT.*mCaT.^3.*hCaT -gCaS.*mCaS.^3.*hCaS -gA.*mA.^3.*hA -gH.*mH -gKCa.*mKCa.^4 -gKd.*mKd.^4 -gleak -gsyn.*[msyn(2) msyn(1)]);
    sum_Eg = -(-gapp(tt,:).*Vsyn-gNa.*mNa.^3.*hNa.*VNa -gCaT.*mCaT.^3.*hCaT.*VCa -gCaS.*mCaS.^3.*hCaS.*VCa -gA.*mA.^3.*hA.*VK -gH.*mH.*Vh -gKCa.*mKCa.^4.*VK -gKd.*mKd.^4.*VK -gleak.*Vleak -gsyn.*[msyn(2) msyn(1)].*Vsyn_HCO);
	 
    ICa = (gCaT.*mCaT.^3.*hCaT + gCaS.*mCaS.^3.*hCaS) .*(V(tt,:)-VCa);
    tau_V = C./sum_g;  %Membrane time constant.
	 
	%V_inf is steady-state voltage.
	%V_inf = (sum_Eg + Iapp(tt,:))./sum_g;
    V_inf = sum_Eg./sum_g;
	
	%evolve the voltage using exponential Euler
	V(tt+1,:) = V_inf + (V(tt,:)-V_inf).*exp(-dt./tau_V);
end
%%% Figures
V = V(1:end-1,:);

if plot_traces
    figure
    subplot(2,1,1);
    plot(tot_T,V(:,1),'k');xlim([2000 3000])
    subplot(2,1,2);plot(tot_T,V(:,2),'r');
    xlim([2000 3000])
end