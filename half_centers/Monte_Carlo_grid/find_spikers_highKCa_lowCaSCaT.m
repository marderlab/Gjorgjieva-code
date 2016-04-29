
spikers = find(bp(pairs(record2),1)<2 & bp(pairs(record2,2),1)<2);
vecc = pairs(record2(spikers),:);

clear gg_KCa gg_CaT gg_CaS
for ii=1:size(spikers,1)
    for jj=1:2
        gg_KCa(ii,jj) = g(vecc(ii,jj),6);
        gg_CaT(ii,jj) = g(vecc(ii,jj),3);
        gg_CaS(ii,jj) = g(vecc(ii,jj),4);
    end
end

collect = [];
collect_sp = [];
for jj=1:size(spikers,1)
    if ((gg_KCa(jj,:)>[100 100]) & (gg_CaT(jj,:)<[1.8 1.8]) & (gg_CaS(jj,:)<[1.8 1.8]))
        collect = [collect jj];
        collect_sp = [collect_sp spikers(jj)];
    end
end

%to get the pairs:
vecc(collect,:)

pairs(record2(collect_sp),:)