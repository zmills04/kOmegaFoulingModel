nX = 64;
nY = 100;
Re = 10;
visc = 0.25;
alpha = visc/0.71;

Umax = 1.5*Re*visc/nY*2;
Yvals = (0.5:(nY-0.5)) - nY/2;
Uvals = Umax * (1. - Yvals.*Yvals/(nY/2)/(nY/2));

Uvals = [0,0,Uvals,0,0];
nY += 4;
Tin = 2;
A = zeros(nX*nY);
T = zeros(nX*nY,1);
b = zeros(nX*nY,1);
for j = 3:nY-2
  for i = 1:nX
    ind = i + (j-1)*nX;
    A(ind,ind) = 4*alpha;
       
    if(i > 1)
      A(ind,ind-1) = -(Uvals(j)/2 + alpha);
    else
      b(ind) = Tin*(Uvals(j)/2 + alpha);
    end
    
    if(i < nX)
      A(ind,ind+1) = (Uvals(j)/2 - alpha);
    else
      A(ind,ind) += (Uvals(j)/2 - alpha);
    end
    
    if(j > 3)
      A(ind,i+nX*(j-2)) = -alpha;
    end

    if(j < nY-2)
      A(ind,i+nX*(j)) = -alpha;
    end
  end
end


Tnew = pcg(A,b);  
Tnew2 = reshape(Tnew,nX,nY);

















