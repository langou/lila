%
   function [ A ] = lapack_geqr2( m, n, i, A )
%
   ilo = i;
   ihi = ilo + n - 1;
%
   for j = ilo:ihi, 
%B = A(ilo:m,ilo:ihi)
      [ A(j:m,j) ] = lapack_larfg( A(j:m,j) );
%[B-A(ilo:m,ilo:ihi),ones(m-ilo+1,1),A(ilo:m,ilo:ihi)]
      [ A(j:m,j+1:ihi) ] = lapack_larfL( A(j:m,j), A(j:m,j+1:ihi) );
%C = A(ilo:m,ilo:ihi)
   end
%
   end
%
