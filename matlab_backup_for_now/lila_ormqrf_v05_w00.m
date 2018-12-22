   function [ A ] = lila_ormqrf_v05_w00( m, n, k, i, j, mt, A, T )
%
   for ii = i:i+k-1, 
%
   V = tril(A(ii:m,ii), -1) + eye(m-ii+1,1);
%
   H = (eye(m-ii+1,m-ii+1) - V * ( lapack_larft( V ) * V' ) );
%
   A(ii:m,j:j+n-1) = H*A(ii:m,j:j+n-1);
%
   end
%
