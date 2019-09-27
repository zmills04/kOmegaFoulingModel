loc = Par1.loc;
[~,lt] = sort(loc);

Par1.pos(:,1) = Par1.pos(lt,1);
Par1.pos(:,2) = Par1.pos(lt,2);
Par1.loc = Par1.loc(lt);
Par1.timer = Par1.timer(lt);
Par1.depTimer = Par1.depTimer(lt);
Par1.numRep = Par1.numRep(lt);
Par1.type = Par1.type(lt);









% loc = load('Par_Loc.txt');
inds = 0:length(loc)-1;

[sinds,lt] = sort(loc);
indsSorted = inds(lt);
n = 256;
locNew = zeros(n,1);

for i = 1:n
    locNew(i) = loc(indsSorted(i)+1);
end

sinds = [sinds(1:n),locNew];



return

load SortLocs.txt

% startVals = SortLocs(1:256:length(SortLocs));


indsSort = zeros(length(loc),1);

nSorts = length(loc)/256;

for i = 0:nSorts-1
    locTemp = loc(((i*256)+1):((i+1)*256));
    [~,lt] = sort(locTemp);
    lt = lt + i*256;
    indsSort(((i*256)+1):((i+1)*256)) = inds(lt);
end

indsSort = [indsSort,SortLocs];






