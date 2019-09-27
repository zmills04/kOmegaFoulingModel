RevDir = [2, 1, 4, 3, 6, 5, 8, 7];
Cx = [1,-1,0,0,1,-1,1,-1];
Cy = [0,0,1,-1,1,-1,-1,1];
[x,y] = size(ibbssort);
ibbssort_errs = zeros(x,y+6);
errind = 1;
for i = 1:length(ibbssort)
    xval = ibbssort(i,1);
    yval = ibbssort(i,2);
    dir = ibbssort(i,3);
    
    xrev = xval;
    yrev = yval;
    dirrev = RevDir(dir);
    
    xrevneigh = xval + Cx(dirrev);
    yrevneigh = yval + Cy(dirrev);
    dirrevneigh = dirrev;
    corr_vals = ones(1,6)*-1;
    err = 0;
    if(xrev ~= ibbssort(i,4))
        err = 1;
        corr_vals(1) = xrev;
    end
    if(yrev ~= ibbssort(i,5))
        err = 1;
        corr_vals(2) = yrev;
    end
    if(dirrev ~= ibbssort(i,6))
        err = 1;
        corr_vals(3) = dirrev;
    end
    if(xrevneigh ~= ibbssort(i,7))
        err = 1;
        corr_vals(4) = xrevneigh;
    end
    if(yrevneigh ~= ibbssort(i,8))
        err = 1;
        corr_vals(5) = yrevneigh;
    end
    if(dirrevneigh ~= ibbssort(i,9))
       err = 1;
       corr_vals(6) = dirrevneigh;
    end
    if(err)
        ibbssort_errs(errind,:) = [ibbssort(i,:),corr_vals];
        errind = errind+1;
    end
end
