function [ x, arnoldi_res, true_res, orth_level, repres_level, cond_level ] = gmres__( A, x, b, Ml, Mr, m, max_it, tol )
%
   hh_lvl1 = 1;
   hh_lvl2 = 0;
   mgs_lvl1 = 0;
%
   if ( hh_lvl1 + hh_lvl2 + mgs_lvl1 ~= 1 ) fprintf('ooops\n'); return; end;
%
   check = 1;
%
   [n] = size(A,1);
%
   H(1:m+1,1:m) = zeros(m+1,m);
   cs(1:m) = zeros(m,1);
   sn(1:m) = zeros(m,1);
   s = zeros(m+1,1);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (mgs_lvl1), 
   Q(1:n,1:m+1) = zeros(n,m+1);
   if ( check ), H_(1:m+1,1:m) = zeros(m+1,m); end;
end
%
if (hh_lvl1), 
   tau = zeros(m,1);
   V(1:n,1:m+1) = zeros(n,m+1);
   if ( check ), Q(1:n,1:m+1) = zeros(n,m+1); end;
   if ( check ), H_(1:m+1,1:m) = zeros(m+1,m); end;
end
%
if (hh_lvl2), 
   T = zeros(m,m);
   V(1:n,1:m+1) = zeros(n,m+1);
   if ( check ), Q(1:n,1:m+1) = zeros(n,m+1); end;
   if ( check ), H_(1:m+1,1:m) = zeros(m+1,m); end;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
   iter = 1;
%
   bnrm2 = norm(Ml*b,2);
%
   r = Ml * ( b-A*Mr*x);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (mgs_lvl1), 
   Q(:,1) = r;
   [ Q(:,1), s(1) ] = orth_mgs_lvl1( Q(:,1) );
end
%
if (hh_lvl1), 
   V(:,1) = r;
   [ V(:,1), s(1), tau ] = orth_hh_lvl1( V(:,1), tau );
end
%
if (hh_lvl2), 
   Q(:,1) = r;
   [ V(:,1), s(1), T(1,1), Q(:,1) ] = orth_hh_lvl2( [], Q(:,1), [] );
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   arnoldi_res(1) = abs( s(1) ) / bnrm2;
   true_res(1) = abs( s(1) ) / bnrm2;
true_res(1)
%  if ( arnoldi_res(1) < tol ), x=Mr*x; return; end
%
%  checks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (hh_lvl1), 
   [ Q(:,1) ] = orth_hh_lvl1_giveq( V(:,1), tau );
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   beta = s(1);
   orth_level(iter) = norm( eye(1) - Q(1:n,1)'*Q(1:n,1), 'fro' );
   repres_level(iter) = norm( r - Q(1:n,1:1)*beta, 'fro' ) / norm( beta, 'fro' );
   cond_level(iter) = cond( [ r ] );

   fprintf('iteration = %3d, iteration = (%3d,%3d), orth level = %6.2e, repres level = %6.2e, cond level = %6.2e, arnoldi res = %6.2e, true res = %6.2e\n', iter-1, 0, 0, orth_level(iter), repres_level(iter), cond_level(iter), arnoldi_res(iter), true_res(iter));
%
   j = 0;
   while (( iter <= max_it+1 )&&(arnoldi_res(iter)>tol)),

      i = 0;
      j = j+1;
      while (( i < m )&( iter <= max_it+1 )&&(arnoldi_res(iter)>tol)),

         i = i+1;
         iter = iter+1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (mgs_lvl1), 
         Q(:,i+1) = Ml*(A*Mr)*Q(:,i);
         [ Q(:,1:i+1), H(1:i+1,i) ] = orth_mgs_lvl1( Q(:,1:i+1) );
end
%
if (hh_lvl1), 
         [ V(:,i+1) ] = orth_hh_lvl1_giveq( V(:,1:i), tau );
         V(:,i+1) = Ml*(A*Mr)*V(:,i+1);
         [ V(:,1:i+1), H(1:i+1,i), tau ] = orth_hh_lvl1( V(:,1:i+1), tau );
end
%
if (hh_lvl2), 
         V(:,i+1) = Q(:,i);
         V(:,i+1) = Ml*(A*Mr)*V(:,i+1);
         [ V(:,i+1), H(1:i+1,i), T(1:i+1,i+1), Q(:,i+1) ] = orth_hh_lvl2( V(:,1:i), V(:,i+1), T(1:i,1:i) );
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        checks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (hh_lvl1), 
         [ Q(:,i+1) ] = orth_hh_lvl1_giveq( V(:,1:i+1), tau );
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         H_(1:i+1,i) = H(1:i+1,i);
         orth_level(iter) = norm( eye(i+1) - Q(1:n,1:i+1)'*Q(1:n,1:i+1), 'fro' );
         repres_level(iter) = norm( [ r Ml*(A*Mr)*Q(1:n,1:i) ] - Q(1:n,1:i+1)*[ [beta; zeros(i,1)], H_(1:i+1,1:i) ], 'fro' ) / norm( [[beta; zeros(i,1)], H_(1:i+1,1:i) ], 'fro' );
         cond_level(iter) = cond( [ r Ml*(A*Mr)*Q(1:n,1:i) ] );

         for k = 1:i-1,
            temp     =  cs(k)*H(k,i) + sn(k)*H(k+1,i);
            H(k+1,i) = -sn(k)*H(k,i) + cs(k)*H(k+1,i);
            H(k,i)   = temp;
         end
         if ( H(i+1,i) == 0.0 ),
            cs(i) = 1.0;
            sn(i) = 0.0;
         elseif ( abs(H(i+1,i)) > abs(H(i,i)) ),
            temp = H(i,i) / H(i+1,i);
            sn(i) = 1.0 / sqrt( 1.0 + temp^2 );
            c = temp * sn(i);
         else
            temp = H(i+1,i) / H(i,i);
            cs(i) = 1.0 / sqrt( 1.0 + temp^2 );
            sn(i) = temp * cs(i);
         end
         temp   = cs(i)*s(i);
         s(i+1) = -sn(i)*s(i);
         s(i)   = temp;
         H(i,i) = cs(i)*H(i,i) + sn(i)*H(i+1,i);
         H(i+1,i) = 0.0;
         arnoldi_res(iter)  = abs(s(i+1)) / bnrm2;

%        checks
         y__ = H(1:i,1:i) \ s(1:i);
         x__ = x + Q(:,1:i)*y__;
         r__ = Ml*( b-(A*Mr)*x__ );
         true_res(iter) = norm( r__, 2 ) / bnrm2;
         fprintf('iteration = %3d, iteration = (%3d,%3d), orth level = %6.2e, repres level = %6.2e, cond level = %6.2e, arnoldi res = %6.2e, true res = %6.2e\n', iter-1, j, i, orth_level(iter), repres_level(iter), cond_level(iter), arnoldi_res(iter), true_res(iter));

      end

      y = H(1:i,1:i) \ s(1:i);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (mgs_lvl1), 
     x = x + Q(:,1:i) * y(1:i);
end
%
if (hh_lvl1), 
      x1 = [y(1:i);zeros(n-i,1)];
      alpha = tau(i) * x1(i,1);
      x1(i,1) = x1(i,1) - V(i,i) * alpha;
      x1(i+1:n,1) = - V(i+1:n,i) * alpha;
      for k = i-1:-1:1,
         alpha = tau(k) * ( V(k:n,k)' * x1(k:n,1) );
         x1(k:n,1) = x1(k:n,1) - V(k:n,k) * alpha;
      end 
      x = x + x1;
end
%
if (hh_lvl2), 
     x = x + Q(:,1:i) * y(1:i);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


      if (( iter <= max_it+1 )&&(arnoldi_res(iter)>tol)),
         r = Ml*( b-(A*Mr)*x );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (mgs_lvl1), 
         Q(:,1) = r;
         [ Q(:,1), s(1) ] = orth_mgs_lvl1( Q(:,1) );
end
%
if (hh_lvl1), 
         V(:,1) = r;
         [ V(:,1), s(1), tau ] = orth_hh_lvl1( V(:,1), tau );
end
%
if (hh_lvl2), 
         V(:,1) = r;
         [ V(:,1), s(1), T(1,1), Q(:,1) ] = orth_hh_lvl2( [], V(:,1), [] );
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%  checks
if (hh_lvl1), 
   [ Q(:,1) ] = orth_hh_lvl1_giveq( V(:,1), tau );
end
         beta = s(1);

      end

   end