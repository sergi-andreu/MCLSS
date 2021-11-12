"""
    Q,H=arnoldi(A,b,m)

A simple implementation of the Arnoldi method.
The algorithm will return an Arnoldi "factorization":
Q*H[1:m+1,1:m]-A*Q[:,1:m]=0
where Q is an orthogonal basis of the Krylov subspace
and H a Hessenberg matrix.

The function `my_hw1_gs(Q,w,k)` needs to be available.

Example:
```julia-repl
using Random
A=randn(100,100); b=randn(100);
m=10;
Q,H=arnoldi(A,b,m);
println("1:st should be zero = ", norm(Q*H-A*Q[:,1:m]));
println("2:nd Should be zero = ", norm(Q'*Q-I));
```

"""
function arnoldi(A,b,m)

    n=length(b);
    Q=zeros(n,m+1);
    H=zeros(m+1,m);
    Q[:,1]=b/norm(b);

    for k=1:m
        w=A*Q[:,k]; # Matrix-vector product with last element
        # Orthogonalize w against columns of Q.
        # Implement this function or replace call with code for orthogonalizatio
        h,β,z=my_hw1_gs(Q,w,k);
        #Put Gram-Schmidt coefficients into H
        H[1:(k+1),k]=[h;β];
        # normalize
        Q[:,k+1]=z/β;
    end
    return Q,H
end


function CGS(Q, b, k)
    
    Q = Q[:,1:k];
    h = Q'*b;
    z = b - Q*h;
    beta = norm(z);
    return h, beta, z
end
    
function MGS(Q, b, k)
    
    z = b; h = zeros(k,1);
    
    for i=1:k
        h[i] = Q[:,i]'*b;
        z = z - h[i]*Q[:,i];
    end
    
    beta = norm(z)
    return h, beta, z
end

function DGS(Q, b, k)
    Q = Q[:,1:k]
    
    h = Q'*b;
    z = b - Q*h;
    
    g = Q'*z;
    z = z - Q*g;
    
    h = h + g;
    beta = norm(z);
    
    return h, beta, z
end


function TGS(Q, b, k)
    Q = Q[:,1:k];
    
    h = Q'*b;
    z = b-Q*h;
    
    g = Q'*z;
    z = z - Q*g;
    
    y = Q'*z;
    z = z- Q*y;
    
    h = h+g+y;
    beta = norm(z)
    
    return h, beta, z
end
    

