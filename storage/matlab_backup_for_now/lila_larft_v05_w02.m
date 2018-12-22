   function [ T ] = lila_larft_v05_w02( m, n, i, mt, A, T )
%
   T( i:i+n-1,i:i+n-1 ) = lapack_larft( A( i:m,i:i+n-1 ) );
%
   end
