%
   function [ Q ] = lapack_orgqr2( m, n, i, A, Q )
%
      Q(i:m,i:i+n-1) = eye( size ( Q(i:m,i:i+n-1) ) );
      for j = i+n-1:-1:i,
         [ Q(j:m,j:i+n-1) ] = lapack_larfL( A(j:m,j), Q(j:m,j:i+n-1) );
      end
%
   end

