figure;set(gcf,'position',[0 0 2000 400]);

result = [];
for ii = 1:40
    ii
    for gr=1:10
        subplot(2,5,gr);plot(tot_T,VV{ii,gr}(:,1),'k');
        hold on;plot(tot_T,VV{ii,gr}(:,2),'r');hold off;
        xlim([1500 2000])
    end
    prompt = 'Enter the number of sim?';
    x = input(prompt);
    result = [result x];
end
