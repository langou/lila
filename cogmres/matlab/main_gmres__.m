%
   clear
   n = 31;
%
   A = zeros(n,1);
   for i = 1:n, A(i,i) = 1.0e+1; end
   for i = 1:n-1, A(i,i+1) = 1.0e+00; end
   for i = 1:n-1, A(i+1,i) = 2.0e+00; end
   b = (n/2)*ones(n,1)-(1:n)';
   Ml = eye(n)+1e-1*diag(1:n);
Ml=eye(n);
   Mr = 1e-1*diag(1:n);
   %A = A*Mr;
   %Mr = eye(n);
   restrt = 10;
   max_it = 15;
   tol = 1e-10;
   x = zeros(n,1);
%
   [ x, arnoldi_res, true_res, orth_level, repres_level, cond_level ] ...
         = gmres__( A, x, b, Ml, Mr, restrt, max_it, tol );
%
   fprintf('                                                                                                                                          check = %6.2e\n', norm(Ml*( b-A*Mr*x ),2)/norm(Ml*b,2));

