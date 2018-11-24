%
   clear
   n = 10;
%
   A = zeros(n,1);
   for i = 1:n, A(i,i) = 1.0e+1; end
   for i = 1:n-1, A(i,i+1) = 1.0e+00; end
   for i = 1:n-1, A(i+1,i) = 2.0e+00; end
   x = ones(n,1)-1e-1*(1:n)';
   M = eye(n)+1e-1*diag(1:n);
   restrt = 3;
   max_it = 10;
   tol = 1e-4;
   b = ones(n,1);
%
   [ x, arnoldi_res, true_res, orth_level, repres_level ] ...
         = gmres__( A, x, b, M, restrt, max_it, tol );
%
   fprintf('                                                                                                                   check = %6.2e\n', norm(M\( b-A*x ),2)/norm(b,2));

