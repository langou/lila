%
   function [ a ] = lapack_larfg( a )
%
      m = size(a,1);
      norma = norm( a, 2);
      if ( a(1) > 0 ) a(1) = a(1) + norma; else a(1) = a(1) - norma; end
      a(2:m,1) = a(2:m) / a(1);
%
      if ( a(1) > 0 ) a(1) = - norma; else a(1) = norma; end
%
   end
