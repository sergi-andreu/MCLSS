function H=amazing_hessenberg_red(A);
% A more refined way of carrying out 
% the hessenberg reduction 
% (Algorithm 2 from Block 3)

n=length(A);

for k=1:n-2
    x=A(k+1:end,k);
    ei = zeros(length(x),1); ei(1) = 1;
    rho = 1; alpha = rho*norm(x);
    z=x-alpha*ei;
    u=z/norm(z);
    
    A(k+1:n,k:n) = A(k+1:n,k:n) - 2*u*(u'*A(k+1:n,k:n)); 
    A(1:n,k+1:n) = A(1:n,k+1:n) - 2*(A(1:n, k+1:n)*u)*u';
    
end
H=A; % should be a hessenberg matrix with same eigenvalues as input A