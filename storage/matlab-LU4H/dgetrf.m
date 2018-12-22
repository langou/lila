%
   clear
%
%  Hopefully named the same in LAPACK 
%
%
   m = 21;
   n = 7;
   A = randn(m,n);
   As = A;
   nrmA = norm(A, 'fro');
%
   for k = 1:n,
      A(k+1:m,k) = A(k+1:m,k) / A(k,k);
      A(k+1:m,k+1:n) = A(k+1:m,k+1:n) -  A(k+1:m,k) * A(k,k+1:n);
   end
%
   U = triu(A); U=U(1:n,1:n);
   L = tril(A,-1)+eye(m,n);
   fprintf('|| A - L*U || / ||A|| = %d', norm(As - L*U,'fro') / nrmA );
   fprintf('\n');




