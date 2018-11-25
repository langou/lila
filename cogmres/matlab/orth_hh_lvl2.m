%
   function [ v, r, t, a ] = orth_hh_lvl2( V, a, T )
%
%  this is Householder using Level 2 BLAS
%  this is the "standard" code that would use 4 synchronizations
%
      m = size(a,1);
      j = size(V,2)+1;
%
      v = zeros(m,1);
      r = zeros(j,1);
      t = zeros(j,1);
      tmp = zeros(j,1);
%
      if (j > 1 )
      tmp(1:j-1,1) = V(1:m,1:j-1)' * a(1:m,1);   %%%%%%%%% <- 1st synchronization
      tmp(1:j-1,1) = T(1:j-1,1:j-1)' * tmp(1:j-1,1) ;
      a(1:m,1) = a(1:m,1) - V(1:m,1:j-1) * tmp(1:j-1,1) ;
%
      r(1:j-1,1) = a(1:j-1,1);                   %%%%%%%%% <- there is a broadcast here if we want all procs to have the Hessenberg matrix
      end
%
      normx = norm( a(j+1:m,1) , 2);             %%%%%%%%% <- 2nd synchronization
      norma = sqrt( normx*normx +  a(j,1)*a(j,1) );
      if ( a(j,1) > 0.0e+00 ) v(j,1) = a(j,1) + norma; else v(j,1) = a(j,1) - norma; end
      v(j+1:m,1) = a(j+1:m,1) / v(j,1);
      if ( a(j,1) > 0.0e+00 ) r(j,1) = - norma; else r(j,1) = norma; end
      tau(j) = 2.0e+00 / ( 1.0e+00 + ( normx / v(j,1) )^2 );
      v(j,1) = 1.0e+00;
%
      if (j > 1 )
      t(1:j-1,1) = V(1:m,1:j-1)'*v(1:m,1);       %%%%%%%%% <- 3rd synchronization
      t(1:j-1,1) = -t(1:j-1,1)*tau(j);
      t(1:j-1,1) = T(1:j-1,1:j-1)*t(1:j-1,1);
      end
      t(j,1) = tau(j);
%
      a(j,1) = - tau(j) ;
      a(j+1:m,1) = v(j+1:m,1) * a(j,1) ;
      a(j,1) = 1.0e+00 + a(j,1) ;
      if (j > 1 )
      a(1:j-1,1) = 0.e+00;
      tmp(1:j-1,1) = V(1:m,1:j-1)' * a(1:m,1) ;  %%%%%%%%% <- 4th synchronization
      tmp(1:j-1,1) = T(1:j-1,1:j-1) * tmp(1:j-1,1) ;
      a(1:m,1) = a(1:m,1) - V(1:m,1:j-1) * tmp(1:j-1,1) ;
      end
%
   end
%
