RevDir = [1, 0, 3, 2, 5, 4, 7, 6];
DIST_SIZE = 417792;
CHANNEL_LENGTH_FULL = 2^nextpow2(2000);
ibbNeigh = [DIST_SIZE - 1, -DIST_SIZE + 1,...
    DIST_SIZE - CHANNEL_LENGTH_FULL, -DIST_SIZE + CHANNEL_LENGTH_FULL,...
    DIST_SIZE - (CHANNEL_LENGTH_FULL + 1), -DIST_SIZE + (CHANNEL_LENGTH_FULL + 1),...
	DIST_SIZE - CHANNEL_LENGTH_FULL + 1, -DIST_SIZE + CHANNEL_LENGTH_FULL - 1];

load ibbArr.txt;
zeroInds = find(ibbArr == 0);
ibbArr(zeroInds) = [];

dir = floor(ibbArr/DIST_SIZE);
linLoc = mod(ibbArr,DIST_SIZE);
yLoc = floor(linLoc/CHANNEL_LENGTH_FULL);
xLoc = mod(linLoc,CHANNEL_LENGTH_FULL);
Ibb = [xLoc,yLoc,dir];
[IbbSort,indSort] = sortrows(Ibb,[1,2,3]);

revDirArr = zeros(length(ibbArr),1);
revNeighArr = revDirArr;
for i = 1:length(ibbArr)
    revDirArr(i) = RevDir(dir(i));
    revNeighArr(i) = ibbNeigh(dir(i));
end
    

revDir = linLoc + (revDirArr + 1) * DIST_SIZE;
revNeigh = ibbArr + revNeighArr;

dir1 = floor(revDir/DIST_SIZE);
linLoc1 = mod(revDir,DIST_SIZE);
yLoc1 = floor(linLoc1/CHANNEL_LENGTH_FULL);
xLoc1 = mod(linLoc1,CHANNEL_LENGTH_FULL);

Ibb1 = [xLoc1,yLoc1,dir1];
IbbSort1 = Ibb2(indSort,:);

dir2 = floor(revNeigh/DIST_SIZE);
linLoc2 = mod(revNeigh,DIST_SIZE);
yLoc2 = floor(linLoc2/CHANNEL_LENGTH_FULL);
xLoc2 = mod(linLoc2,CHANNEL_LENGTH_FULL);

Ibb2 = [xLoc2,yLoc2,dir2];
IbbSort2 = Ibb2(indSort,:);

Ibb3 = [IbbSort,IbbSort1,IbbSort2];

% Test inds
Cx = [1,-1,0,0,1,-1,1,-1];
Cy = [0,0,1,-1,1,-1,-1,1];
for i = 1:length(ibbArr)
    xval = Ibb3(i,1);
    yval = Ibb3(i,1);
    dir = Ibb3(i,1);
    
    xrev = xval;
    yrev = yval;
    dirrev = RevDir(dir);
    
    xrevneigh = xval + Cx(dirrev+1);
    yrevneigh = yval + Cy(dirrev+1);
    dirrevneigh = dirrev;
    
    if(xrev ~= Ibb3(i,4))
        fprintf("error in index %d, row %d",i,4);
    end
    if(yrev ~= Ibb3(i,5))
        fprintf("error in index %d, row %d",i,5);
    end
    if(dirrev ~= Ibb3(i,6))
        fprintf("error in index %d, row %d",i,6);
    end
    if(xrevneigh ~= Ibb3(i,7))
        fprintf("error in index %d, row %d",i,7);
    end
    if(yrevneigh ~= Ibb3(i,8))
        fprintf("error in index %d, row %d",i,8);
    end
    if(dirrevneigh ~= Ibb3(i,9))
        fprintf("error in index %d, row %d",i,9);
    end

end





