Acsr = read_sparse_triplets("CSR",64*104);
T = ones(104*64,1);
Bvec = load("Temp_bvec.txt");
Afull = full(Acsr);
%[l,u,p] = ilu(A,struct("droptol",1e-3));
Tcsr = pcg(Acsr,Bvec);


%Bvec = reshape(Bvec,512,104);
%Tnew = A*T;
Tcsr2 = reshape(Tcsr,64,104);
imagesc(Tcsr2)

