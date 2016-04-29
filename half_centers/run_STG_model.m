function run_STG_model
%%% Parameters
C = 0.628; %muF/cm2
VNa = 50; VK = -80; VCa = 80; Vleak = -50; Vsyn = -80;%mV
gleak = 0.01; %mS/cm2
gNa = 700; gCaT = 2; gCaS = 4; 
gKd = 70; gKCa = 40; gA = 50; 
gsyn = 0; tausyn = 20;

Iapp = 0;
Iappstep = 0;

%%% Initial values
V(1,1:2) = [-50 -50];%-70;
Ca(1,1:2) = 0.5*[1 1];%0.2;
mNa = Xinf(V(1,:),25.5,-5.29);
hNa = Xinf(V(1,:),48.9,5.18);
mCaT = Xinf(V(1,:),27.1,-7.2);
hCaT = Xinf(V(1,:),32.1,5.5);
mCaS = Xinf(V(1,:),33,-8.1);
hCaS = Xinf(V(1,:),60,6.2);
mA = Xinf(V(1,:),27.2,-8.7);
h = Xinf(V(1,:),56.9,4.9);
mKCa = m_KCainf(V(1,:),Ca(1,:));
mKd = Xinf(V(1,:),12.3,-11.8);
msyn = Xinf(V(1,:),45,-2);

%%% Simulation
tinit = 0;
tfinal = 1000;
dt = 0.01;
total_time = tinit:dt:tfinal;

for tt=1:length(total_time)
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
    msyninf = Xinf(V(tt,:),45,-2);
    
    I_ion = (1./C) * (-gNa.*mNa.^3.*hNa.*(V(tt,:)-VNa) -gCaT.*mCaT.^3.*hCaT.*(V(tt,:)-VCa) -gCaS.*mCaS.^3.*hCaS.*(V(tt,:)-VCa) -gA.*mA.^3.*hA.*(V(tt,:)-VK) -gKCa.*mKCa.^4.*(V(tt,:)-VK) -gKd.*mKd.^4.*(V-VK) -gleak.*(V(tt,:)-Vleak) -gsyn.*[msyn(2); msyn(1)].*(V(tt,:)-Vsyn) + Iapp + (Iappstep*(t>2000 && t<4000)));
    V(tt+1,:) = V(tt,:) + dt/C*(-I_ion);
    
    mNa = mNa + (1-exp(-dt/taumNa))*(mNainf - mNa);
    hNa = hNa + (1-exp(-dt/tauhNa))*(hNainf - hNa);
    mCaT = mCaT + (1-exp(-dt/taumCaT))*(mCaTinf - mCaT);
    hCaT = hCaT + (1-exp(-dt/tauhCaT))*(hCaTinf - hCaT);
    mCaS = mCaS + (1-exp(-dt/taumCaS))*(mCaSinf - mCaS);
    hCaS = hCaS + (1-exp(-dt/tauhCaS))*(hCaSinf - hCaS);
    mA = mA + (1-exp(-dt/taumA))*(mAinf - mA);
    hA = hA + (1-exp(-dt/tauhA))*(hAinf - hA);
    mKCa = mKCa + (1-exp(-dt/taumKCa))*(mKCainf - mKCa);
    mKd = mKd + (1-exp(-dt/taumKd))*(mKdinf - mKd);
    msyn = msyn + (1-exp(-dt/tausyn))*(msyninf - msyn);

%%% Figures
figure
plot(total_time,V(:,1),'k')
hold on;plot(total_time,V(:,2),'r--');
axis([0 tfinal -100 60])
hold off