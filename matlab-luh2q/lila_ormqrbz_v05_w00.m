%
   function [ Q ] = lila_ormqrbz_v05_w00( m, n, k, i, j, mt, A, T, Q )
%
      Q(i:i+k-1,j:j+n-1) = zeros( size( Q(i:i+k-1,j:j+n-1) ));
%
      for ii = i+k-1:-1:i,
%
         V = tril(A(ii:m,ii), -1) + eye(size(A(ii:m,1)));
%
         H = (eye(m-ii+1,m-ii+1) - V * ( lapack_larft( V ) * V' ) );
%
         Q(ii:m,j:j+n-1) = H * Q(ii:m,j:j+n-1);
%
      end
%
   end
