function [ x, arnoldi_res, true_res, orth_level, repres_level ] = gmres__( A, x, b, M, m, max_it, tol )

   iter = 1;

   bnrm2 = norm(b,2);

   r = M \ ( b-A*x );
   beta = norm(r,2);
   arnoldi_res(1) = beta / bnrm2;
   true_res(1) = beta / bnrm2;
   if ( arnoldi_res(1) < tol ) return, end

   [n,~] = size(A);
   V(1:n,1:m+1) = zeros(n,m+1);
   H(1:m+1,1:m) = zeros(m+1,m);
   cs(1:m) = zeros(m,1);
   sn(1:m) = zeros(m,1);
   s = zeros(m+1,1);

   beta = norm(r,2);
   V(:,1) = r / beta;
   s(1) = beta;

   orth_level(iter) = norm( eye(1) - V(1:n,1)'*V(1:n,1), 'fro' );
   repres_level(iter) = norm( r - V(1:n,1:1)*beta, 'fro' ) / norm( beta, 'fro' );

   fprintf('iteration = %3d, iteration = (%3d,%3d), orth level = %6.2e, repres level = %6.2e, arnoldi res = %6.2e, true res = %6.2e\n', iter-1, 0, 0, orth_level(iter), repres_level(iter), arnoldi_res(iter), true_res(iter));

   j = 0;
   while (( iter <= max_it+1 )&&(arnoldi_res(iter)>tol)),

      i = 0;
      j = j+1;
      while (( i < m )&( iter <= max_it+1 )&&(arnoldi_res(iter)>tol)),

         i = i+1;
         iter = iter+1;
         w = M\A*V(:,i);

         [ V(:,i+1), H(1:i+1,i) ] = orth_mgs_lvl1( V(:,1:i), w );

%        checks
         H_(1:i+1,i) = H(1:i+1,i);
         orth_level(iter) = norm( eye(i+1) - V(1:n,1:i+1)'*V(1:n,1:i+1), 'fro' );
         repres_level(iter) = norm( [ r M\A*V(1:n,1:i) ] - V(1:n,1:i+1)*[ [beta; zeros(i,1)], H_(1:i+1,1:i) ], 'fro' ) / norm( [[beta; zeros(i,1)], H_(1:i+1,1:i) ], 'fro' );

         for k = 1:i-1,
            temp     =  cs(k)*H(k,i) + sn(k)*H(k+1,i);
            H(k+1,i) = -sn(k)*H(k,i) + cs(k)*H(k+1,i);
            H(k,i)   = temp;
         end
         [cs(i),sn(i)] = rotmat( H(i,i), H(i+1,i) );
         temp   = cs(i)*s(i);
         s(i+1) = -sn(i)*s(i);
         s(i)   = temp;
         H(i,i) = cs(i)*H(i,i) + sn(i)*H(i+1,i);
         H(i+1,i) = 0.0;
         arnoldi_res(iter)  = abs(s(i+1)) / bnrm2;

%        checks
         y__ = H(1:i,1:i) \ s(1:i);
         x__ = x + V(:,1:i)*y__;
         r__ = M\( b-A*x__ );
         true_res(iter) = norm( r__, 2 ) / bnrm2;
         fprintf('iteration = %3d, iteration = (%3d,%3d), orth level = %6.2e, repres level = %6.2e, arnoldi res = %6.2e, true res = %6.2e\n', iter-1, j, i, orth_level(iter), repres_level(iter), arnoldi_res(iter), true_res(iter));

      end

      y = H(1:i,1:i) \ s(1:i);
      beta = s(i+1);

      s = [ zeros(i,1); s(i+1)];
      for k = i:-1:1,
         temp    =  cs(k)*s(k) - sn(k)*s(k+1);
         s(k+1) = sn(k)*s(k) + cs(k)*s(k+1);
         s(k)   = temp;
      end
      r = V(1:n,1:i+1) * s(1:i+1) ;
      x = x + V(:,1:i) * y(1:i);

      V(:,1) = r / beta;
      s(1) = beta;

   end
