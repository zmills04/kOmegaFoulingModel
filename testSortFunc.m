function [loctest, diff] = testSortFunc(locs, inds)

[x,numsort] = size(inds);

loctest = zeros(x,numsort);

for i = 1:numsort
    induse = inds(:,i) + 1;
    loctest(:,i) = locs(induse,i);
end

diff = locs(:,2:end) - loctest;



    
    
    
    
    
    


