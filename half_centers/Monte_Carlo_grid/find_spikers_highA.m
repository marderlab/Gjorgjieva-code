
spikers = find(bp(pairs(record2,1),1)<2 & bp(pairs(record2,2),1)<2);
vecc = pairs(record2(spikers),:);

clear gg_A
for ii=1:size(spikers,1)
    for jj=1:2
        gg_A(ii,jj) = g(vecc(ii,jj),5);
    end
end

collect = [];
collect_sp = [];
for jj=1:size(spikers,1)
    if (gg_A(jj,:)>[127.8 127.8])
        collect = [collect jj];
        collect_sp = [collect_sp spikers(jj)];
    end
end

%to get the pairs:
vecc(collect,:)

pairs(record2(collect_sp),:)