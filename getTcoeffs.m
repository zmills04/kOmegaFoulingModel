tinds = load("temp_Inds.txt");
Amat = load("tempSolver_A.txt");

cC = zeros(2000,200,5);

for i = 1:2000
    for j = 1:200
        jj = j+2;
        ind = i+(jj-1)*2048;
        
        for k = 1:5
            tind_loc = tinds(ind,k);
            if(tind_loc > -1)
                cC(i,j,k) = Amat(tind_loc+1);
            else
                cC(i,j,k) = -1;
            end
        end
    end
end
        
eC = cC(:,:,2);
wC = cC(:,:,3);
nC = cC(:,:,4);
sC = cC(:,:,5);
cC = cC(:,:,1);
