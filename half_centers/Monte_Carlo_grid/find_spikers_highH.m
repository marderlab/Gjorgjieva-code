
spikers = find(bp(pairs(record2),1)<2 & bp(pairs(record2,2),1)<2);
vecc = pairs(record2(spikers),:);

clear g_H
for ii=1:size(spikers,1)
    for jj=1:2
        g_H(ii,jj) = g(vecc(ii,jj),8);
    end
end

collect = [];
collect_sp = [];
for jj=1:size(spikers,1)
    if (g_H(jj,:)>[0.247 0.247])
        collect = [collect jj];
        collect_sp = [collect_sp spikers(jj)];
    end
end

%to get the pairs:
vecc(collect,:)

pairs(record2(collect_sp),:)collect