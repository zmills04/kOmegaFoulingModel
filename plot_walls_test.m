load lsc0.txt
load lsbl.txt
lsbl = lsbl + 1;
close all
hold on
for i = 1:100
    C1 = lsc0(lsbl(i,1),:);
    C2 = lsc0(lsbl(i,2),:);
    
    plot([C1(1),C2(1)], [C1(2),C2(2)],'r')
    tVec = C2-C1;
    nVec = [-tVec(2),tVec(1)];
    plot([C1(1),C1(1)+nVec(1)],[C1(2),C1(2)+nVec(2)],'k')
end