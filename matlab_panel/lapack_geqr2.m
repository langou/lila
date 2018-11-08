%
   function [ A ] = lapack_geqr2( m, n, i, A )
%
      for j = i:i+n-1, 
         [ A(j:m,j) ] = lapack_larfg( A(j:m,j) );
         [ A(j:m,j+1:i+n-1) ] = lapack_larfL( A(j:m,j), A(j:m,j+1:i+n-1) );
      end
%
   end
%
