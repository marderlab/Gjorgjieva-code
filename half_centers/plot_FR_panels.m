

figure
for ii=1:length(rat1)
    for jj=1:length(rat2)
        subplot(3,4,(ii-1)*4+jj);
        plot(mu,squeeze(FR(ii,jj,:,:)),'.-');
        xlim([0 1.5]);ylim([0 100]);
    end
end

figure
for ii=1:length(rat1)
    for jj=1:length(rat2)
        subplot(3,4,(ii-1)*4+jj);
        plot(mu,squeeze(FR_0p1(ii,jj,:,:)),'.-');
        xlim([0 1.5]);ylim([0 100]);
    end
end

figure
for ii=1:length(rat1)
    for jj=1:length(rat2)
        subplot(3,4,(ii-1)*4+jj);
        plot(mu,squeeze(FR_0p2(ii,jj,:,:)),'.-');
        xlim([0 1.5]);ylim([0 100]);
    end
end

figure
for ii=1:length(rat1)
    for jj=1:length(rat2)
        subplot(3,4,(ii-1)*4+jj);
        plot(mu,squeeze(FR_0p3(ii,jj,:,:)),'.-');
        xlim([0 1.5]);ylim([0 100]);
    end
end

figure
for ii=1:length(rat1)
    for jj=1:length(rat2)
        subplot(3,4,(ii-1)*4+jj);
        plot(mu,squeeze(FR_0p4(ii,jj,:,:)),'.-');
        xlim([0 1.5]);ylim([0 100]);
    end
end

figure
for ii=1:length(rat1)
    for jj=1:length(rat2)
        subplot(3,4,(ii-1)*4+jj);
        plot(mu,squeeze(FR_0p5(ii,jj,:,:)),'.-');
        xlim([0 1.5]);ylim([0 100]);
    end
end

%%

figure
for ii=1:length(rat1)
    for jj=1:length(rat2)
        subplot(3,4,(ii-1)*4+jj);
        hold on;
        plot(mu,squeeze(FR(ii,jj,1,:)),'b.-');
        plot(mu,squeeze(FR_0p1(ii,jj,1,:)),'c.-');
        plot(mu,squeeze(FR_0p2(ii,jj,1,:)),'k.-');
        plot(mu,squeeze(FR_0p3(ii,jj,1,:)),'m.-');
        plot(mu,squeeze(FR_0p4(ii,jj,1,:)),'g.-');
        plot(mu,squeeze(FR_0p5(ii,jj,1,:)),'r.-');        
        xlim([0 1.5]);ylim([0 40]);
    end
end

figure
for ii=1:length(rat1)
    for jj=1:length(rat2)
        subplot(3,4,(ii-1)*4+jj);
        hold on;
        plot(mu,squeeze(FR(ii,jj,2,:)),'b.-');
        plot(mu,squeeze(FR_0p1(ii,jj,2,:)),'c.-');
        plot(mu,squeeze(FR_0p2(ii,jj,2,:)),'k.-');
        plot(mu,squeeze(FR_0p3(ii,jj,2,:)),'m.-');
        plot(mu,squeeze(FR_0p4(ii,jj,2,:)),'g.-');
        plot(mu,squeeze(FR_0p5(ii,jj,2,:)),'r.-');        
        xlim([0 1.5]);ylim([0 40]);
    end
end

figure
for ii=1:length(rat1)
    for jj=1:length(rat2)
        subplot(3,4,(ii-1)*4+jj);
        hold on;
        plot(mu,squeeze(FR(ii,jj,3,:)),'b.-');
        plot(mu,squeeze(FR_0p1(ii,jj,3,:)),'c.-');
        plot(mu,squeeze(FR_0p2(ii,jj,3,:)),'k.-');
        plot(mu,squeeze(FR_0p3(ii,jj,3,:)),'m.-');
        plot(mu,squeeze(FR_0p4(ii,jj,3,:)),'g.-');
        plot(mu,squeeze(FR_0p5(ii,jj,3,:)),'r.-');        
        xlim([0 1.5]);ylim([0 50]);
    end
end

figure
for ii=1:length(rat1)
    for jj=1:length(rat2)
        subplot(3,4,(ii-1)*4+jj);
        hold on;
        plot(mu,squeeze(FR(ii,jj,4,:)),'b.-');
        plot(mu,squeeze(FR_0p1(ii,jj,4,:)),'c.-');
        plot(mu,squeeze(FR_0p2(ii,jj,4,:)),'k.-');
        plot(mu,squeeze(FR_0p3(ii,jj,4,:)),'m.-');
        plot(mu,squeeze(FR_0p4(ii,jj,4,:)),'g.-');
        plot(mu,squeeze(FR_0p5(ii,jj,4,:)),'r.-');        
        xlim([0 1.5]);ylim([0 70]);
    end
end

figure
for ii=1:length(rat1)
    for jj=1:length(rat2)
        subplot(3,4,(ii-1)*4+jj);
        hold on;
        plot(mu,squeeze(FR(ii,jj,5,:)),'b.-');
        plot(mu,squeeze(FR_0p1(ii,jj,5,:)),'c.-');
        plot(mu,squeeze(FR_0p2(ii,jj,5,:)),'k.-');
        plot(mu,squeeze(FR_0p3(ii,jj,5,:)),'m.-');
        plot(mu,squeeze(FR_0p4(ii,jj,5,:)),'g.-');
        plot(mu,squeeze(FR_0p5(ii,jj,5,:)),'r.-');        
        xlim([0 1.5]);ylim([0 100]);
    end
end

%%

%diff_FR = zeros(length(rat1),length(rat2),length(rat3));
count = 0;
for ii=1:length(rat1)
    for jj=1:length(rat2)
        for kk=1:length(rat3)
            count = count+1;
            FR_temp(count,:) = FR(ii,jj,kk,:);
            temp{count}=[ii jj kk];
        end
    end
end

count2 = 0;
diff_FR_L1 = zeros(nchoosek(count,2),1);
diff_FR_L2 = zeros(nchoosek(count,2),1);

for ii=1:count
    for jj=ii+1:count
        count2 = count2+1;
        diff_FR_L1(count2) = mean(abs(FR_temp(ii,:)-FR_temp(jj,:)));
        diff_FR_L2(count2) = mean((FR_temp(ii,:)-FR_temp(jj,:)).^2);
        temp2{count2} = [ii jj];
    end
end

figure;plot(diff_FR_L1,'o')
hold on;plot(diff_FR_L2,'r.')