%load FR_mat.mat
mu_vec = [0 0.5 1.0];
sig_vec = [0.0 0.2 0.5];
  
gNa = 80;
rat1 = [0.05 0.1 0.2];
gK = rat1*gNa;
rat2 = [0.05 0.1 0.2 0.3];
gA = rat2*gNa;
rat3 = [0.00001 0.0001 0.001 0.0032 0.01];
gh = rat3*gNa;

FR = FR_0;%FR_0p1;
clear FR_R* FR_D* R_A R_h D_A D_h ind* pair* y*
for rr=1%:length(rat1)
    count = 0;
    for ii=1:length(gA)
        for jj=ii:length(gA)
            count = count+1;
            R_A(count) = gA(ii)/gA(jj); 
            D_A(count) = gA(ii)-gA(jj);         
            pair_A(count,1:2) = [ii,jj];
        end
    end
    count = 0;
    for ii=1:length(gh)
        for jj=ii:length(gh)
            count = count+1;
            R_h(count) = gh(ii)/gh(jj); 
            D_h(count) = gh(ii)-gh(jj);         
            pair_h(count,1:2) = [ii,jj];
        end
    end

    for mm=1:size(pair_A,1)
        for nn=1:size(pair_h,1)
            FR_R(mm,nn) = FR(rr,pair_A(mm,1),pair_h(nn,1),1)/FR(rr,pair_A(mm,2),pair_h(nn,2),1); 
            FR_D(mm,nn) = FR(rr,pair_A(mm,1),pair_h(nn,1),1)-FR(rr,pair_A(mm,2),pair_h(nn,2),1); 
            FR_only(mm,nn) = FR(rr,pair_A(mm,1),pair_h(nn,1),1); 
            FR_only2(mm,nn) = FR(rr,pair_A(mm,2),pair_h(nn,2),1);
        end
    end

    [y_A,ind_A]=sort(R_A);
    [y_h,ind_h]=sort(R_h);
    figure;
    imagesc(abs(FR_D(ind_A,ind_h)))
    set(gca,'ydir','normal')
    set(gca,'xtick',[1:15],'xticklabel',round(y_h*10^3)/10^3)
    set(gca,'ytick',[1:10],'yticklabel',round(y_A*10^3)/10^3)
    colorbar
    set(gca,'fontsize',20)
    title('|\Delta FR (g_A,g_h)|')
    xlabel('g_{h1}/g_{h2}');
    ylabel('g_{A1}/g_{A2}');

    [yd_A,ind_dA]=sort(D_A);
    [yd_h,ind_dh]=sort(D_h);
    figure;
    imagesc(abs(FR_D(ind_dA,ind_dh)));format short
    set(gca,'xtick',[1:15],'xticklabel',round(yd_h*10^3)/10^3)
    set(gca,'ytick',[1:10],'yticklabel',round(yd_A*10^3)/10^3)
    set(gca,'ydir','normal')
    colorbar
    set(gca,'fontsize',20)
    title('|\Delta FR (g_A,g_h)|')
    xlabel('\Delta g_h = g_{h1}-g_{h2}');
    ylabel('\Delta g_A = g_{A1}-g_{A2}');
end

F1=FR_only(ind_dA,ind_dh);
figure;imagesc(F1(7:10,11:15));
colorbar
set(gca,'xtick',[1:5],'xticklabel',round(gh*10^3)/10^3)
set(gca,'ytick',[1:4],'yticklabel',round(gA*10^3)/10^3)
set(gca,'ydir','normal')
set(gca,'fontsize',20)
title('FR (g_A,g_h)')
xlabel('g_h');
ylabel('g_A');
