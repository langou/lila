   function [ T ] = lila_larft_v05_w00( m, n, i, mt, A, T )
%
   for ii=i:i+n-1,
%
   T(1,ii) = lapack_larft( A(ii:m,ii ) );
%
   end
%
