myA = readbin('lbkappaSolver_A',3);
clA = readbin('clA');
myAmat = sparse(myA(:,1),myA(:,2),myA(:,3), 104448, 104448);
clAmat = sparse(myA(:,1),myA(:,2),clA(:), 104448, 104448);
len = 512*204;

myB = getVector('lbkappaSolver_bVec',len);
clB = getVector('clB',len);

myX = getVector('lbkappaSolver',len);
clX = getVector('clX',len);

printf("total diff bVec = %g\n", sum(abs(myB - clB)));
printf("total diff xVec = %g\n", sum(abs(myX - clX)));
printf("total diff Amat = %g\n", sum(abs(myA(:,3) - clA)));

myR = getVector('rVec',len);
clR = getVector('clR',len);
calcR = clB - clAmat*clX;

printf("sum(myR - calcR) = %g, mean(myR - calcR) = %g\n", sum(abs(myR - calcR)), mean(abs(myR - calcR)));
printf("sum(clR - calcR) = %g, mean(clR - calcR) = %g\n", sum(abs(clR - calcR)), mean(abs(clR - calcR)));
printf("sum(myR - clR) = %g, mean(myR - clR) = %g\n", sum(abs(myR - clR)), mean(abs(myR - clR)));

clf
hold on;
plot(myR,'r')
plot(clR,'k')
%plot(calcR,'b')





