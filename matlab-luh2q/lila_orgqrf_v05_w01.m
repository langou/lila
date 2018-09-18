%
   function [ Q ] = lila_orgqrf_v05_w01( m, n, i, mt, A, T, Q )
%
   ilo = i;
   ihi = ilo + n - 1;
%
   Q(ilo:m,ilo:ihi) = eye( size ( Q(ilo:m,ilo:ihi) ) );
   for j = ihi:-1:ilo,
      [ Q(j:m,j:ihi) ] = lapack_larfL( A(j:m,j), Q(j:m,j:ihi) );
   end
%
end

