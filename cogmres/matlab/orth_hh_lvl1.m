%
   function [ V, r, tau, a ] = orth_hh_lvl1( V, tau )
%
      m = size(V,1);
      j = size(V,2);
%
      v = zeros(m,1);
      r = zeros(j,1);
%
      for i = 1:j-1,
         alpha = tau(i) * ( V(i:m,i)' * V(i:m,j) );
         V(i:m,j) = V(i:m,j) - V(i:m,i) * alpha;
      end 
%
      if (j>1), r(1:j-1,1) = V(1:j-1,j); end;
%
      normx = norm( V(j+1:m,j) , 2);
      norma = sqrt( normx*normx +  V(j,j)*V(j,j) );
      if ( V(j,j) > 0.0e+00 ) r(j,1) = - norma; else r(j,1) = norma; end
      if ( V(j,j) > 0.0e+00 ) V(j,j) = V(j,j) + norma; else V(j,j) = V(j,j) - norma; end
      V(j+1:m,j) = V(j+1:m,j) / V(j,j);
      tau(j) = 2.0e+00 / ( 1.0e+00 + ( normx / V(j,j) )^2 );
      V(1:j-1,j) = zeros(j-1,1);
      V(j,j) = 1.0e+00;
%
   end
%
