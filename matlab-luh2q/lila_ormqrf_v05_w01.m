   function [ A ] = lila_ormqrf_v05_w01( m, n, k, i, j, mt, A, T )
%
   Q(i:m,i:i+n-1) = eye( m-i+1, n );
%
   for j = i+n-1:-1:i,
   [ Q(j:m,j:i+n-1) ] = lapack_larfL( A(j:m,i+n-1), Q(j:m,j:i+n-1) );
   end
%
end

