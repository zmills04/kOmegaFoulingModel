function [x1,x2,x3,flags,ress,iters] = testPreconditioner(X, Y, nx, ny, tol, maxit, dt, k, tsteps)

[x1, flag1, res1, iter1] = timeSol(X,Y,nx,ny,k,dt,tsteps,tol,maxit,1);
[x2, flag2, res2, iter2] = timeSol(X,Y,nx,ny,k,dt,tsteps,tol,maxit,2);
[x3, flag3, res3, iter3] = timeSol(X,Y,nx,ny,k,dt,tsteps,tol,maxit,3);
fprintf("For No Preconditioner: total iter = %d, avg res = %g\n", sum(iter1), mean(res1));
fprintf("For Divided: total iter = %d, avg res = %g\n", sum(iter2), mean(res2));
fprintf("For Jacobi Preconditioner: total iter = %d, avg res = %g\n", sum(iter3), mean(res3));
x1 = reshape(x1,nx,ny);
x2 = reshape(x2,nx,ny);
x3 = reshape(x3,nx,ny);

flags = [flag1,flag2,flag3];
ress = [res1,res2,res3];
iters = [iter1,iter2,iter3];
end


function [xvals, flags, ress, iters] = timeSol(X, Y, nx, ny, k, dt, tsteps, tol, maxit, type)

A = zeros(nx*ny,nx*ny);
b = zeros(nx*ny,1);
bConst = b;
xvals = b;
M = A;

dX = X/nx;
dY = Y/ny;
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
        Src = 0;
        
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
        bConst(loc) = Src;
        if(type==3)
            val = A(loc,loc);
            M(loc,loc) = val;
        end
    end
end

if(type==2)
    for i = 1:nx*ny
        val = A(i,i);
        A(i,:) = A(i,:)/val;
    end
    M = diag(A);
end

flags = zeros(tsteps,1);
ress = flags;
iters = flags;

for t = 1:tsteps
    b = bConst + xvals;
    if(type==1)
        [xvals, flags(t), ress(t), iters(t)] = bicg(A, b, tol, maxit);
    elseif(type==3)
        [L,U] = ilu(sparse(A),struct('type','ilutp','droptol',1e-6));
        [xvals, flags(t), ress(t), iters(t)] = bicg(A, b, tol, maxit, L, U);
    else
        b = b./M;
        [xvals, flags(t), ress(t), iters(t)] = bicg(A, b, tol, maxit);
    end
end
end
