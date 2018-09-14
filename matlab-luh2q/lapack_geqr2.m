%
   function [ A ] = lapack_geqr2( m, n, i, A )
%
   ilo = i;
   ihi = ilo + n - 1;
%
   for j = ilo:ihi, 
      [ A(j:m,j) ] = lapack_larfg( A(j:m,j) );
      [ A(j:m,j+1:ihi) ] = lapack_larfL( A(j:m,j), A(j:m,j+1:ihi) );
   end
%
   end
%
