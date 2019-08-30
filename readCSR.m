function [spmat, bVec] = readCSR(fname)

fullname = sprintf("%s_A.txt",fname);
Amat = load(fullname);
fullname = sprintf("%s_IA.txt",fname);
IA = load(fullname);
fullname = sprintf("%s_JA.txt",fname);
JA = load(fullname);
fullname = sprintf("%s_bvec.txt",fname);
bVec = load(fullname);

matsize = length(IA)-1;

row_inds = zeros(length(Amat),1);
col_inds = zeros(length(Amat),1);

[x,y] = size(bVec);
fx = matsize/y;
% bVec = [bVec;zeros(fx-x,y)];
% bVec = reshape(bVec,matsize,1);
bVec = bVec(:,2:end-1);
bVec = reshape(bVec,x*(y-2),1);

cur_ind = 1;
for i = 1:matsize
    numel = IA(i+1) - IA(i);
    if(numel == 0)
        continue;
    end
    
    i_x = mod(i-1,fx);
    i_y = floor((i-1)/fx);
    if(i_x >= x || i_y < 1 || i_y >= 203)
        continue;
    end
    inew = i_x + (i_y-1)*x + 1;
    
    start = IA(i)+1;
    for j = start:IA(i+1)
        j_x = mod(JA(j),fx);
        j_y = floor(JA(j)/fx);
        if(j_x >= x || j_y < 1 || j_y >= 203)
            continue;
        end
        col_inds(cur_ind) = inew;
        jnew = j_x + (j_y-1)*x + 1;
        row_inds(cur_ind) = jnew;
        cur_ind = cur_ind+1;
    end
    
%     start = IA(i)+1;
%     for j = start:IA(i+1)
%         col_inds(cur_ind) = i;
%         row_inds(cur_ind) = JA(j)+1;
%         cur_ind = cur_ind+1;
%     end
end
fprintf("cur_ind = %d, nnz = %d\n",cur_ind,length(Amat));





spmat = sparse(col_inds,row_inds,Amat,length(bVec),length(bVec));



   
    
    









