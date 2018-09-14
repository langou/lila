%
   function [ Q ] = lila_orgqr( m, n, i, A, Q )
%
   ilo = i;
   ihi = ilo + n - 1;
%
   Q(ilo:m,ilo:ihi) = eye( size ( Q(ilo:m,ilo:ihi) ) );
   for j = ihi:-1:ilo,
      [ Q(j:m,j:ihi) ] = larfL( A(j:m,j), Q(j:m,j:ihi) );
   end
%
end

