%
   function [ A ] = lila_geqr2( m, n, i, A )
%
   ilo = i;
   ihi = ilo + n - 1;
%
   for j = ilo:ihi, 
      [ A(j:m,j) ] = larfg( A(j:m,j) );
      [ A(j:m,j+1:ihi) ] = larfL( A(j:m,j), A(j:m,j+1:ihi) );
   end
%
   end
%
