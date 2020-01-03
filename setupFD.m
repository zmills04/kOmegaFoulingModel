function [A,b] = setupFD(X, Y, nx, ny, k, dt, xprev)

A = zeros(nx*ny,nx*ny);
b = zeros(nx*ny,1);
dX = X/nx;
dY = Y/ny;
Tw = 4;
Tn = 3;
Ts = 2;
Te = 1;

Cc = 1+2*k*dt*(1/dX/dX + 1/dY/dY);

for i = 1:nx
    for j = 1:ny
        loc = (i-1)*ny+j;
        A(loc,loc) = Cc;
        Cn = -k*dt/dY/dY;
        Cs = -k*dt/dY/dY;
        Cw = -k*dt/dX/dX;
        Ce = -k*dt/dX/dX;
        Src = xprev(loc);
        
        % adiabatic on west wall
        if(i > 1)
            A(loc, loc-ny) = Cw;
        else
            A(loc,loc) = A(loc,loc) + Cw;
        end
        
        % temp on n,e,s walls
        if(i < nx)
            A(loc, loc+ny) = Ce;
        else
            Src = Src - Ce*Te;
        end
        
        if(j > 1)
            A(loc, loc-1) = Cs;
        else
            Src = Src - Cs*Ts;
        end
        if(j < ny)
            A(loc, loc+1) = Cn;
        else
            Src = Src - Cn*Tn;
        end
        b(loc) = Src;
    end
end

        
        
        

