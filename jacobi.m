function [x] = jacobi(A,b)
n=size(b,1);
x = zeros(n,1);
normVal=Inf; 
%% 
% * _*Tolerence for method*_
tol=1e-5; itr=0;
%% Algorithm: Jacobi Method
%%
while (normVal>tol && itr < 1)
    xold=x;
    sigma = A(1,2:n)*x(2:n);
    x(1) = (1/A(1,1))*(b(1)-sigma);
    for i=2:n-1
        sigma = A(i,[1:i-1,i+1:n])*x([1:i-1,i+1:n]);
        x(i)=(1/A(i,i))*(b(i)-sigma);
    end
    sigma = A(n,1:n-1)*x(1:n-1);
    x(n) = (1/A(n,n))*(b(n)-sigma);
    
    itr=itr+1;
    normVal=abs(xold-x);
end
%%
fprintf('Solution of the system is : \n%f\n%f\n%f\n%f in %d iterations',x,itr);