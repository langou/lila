%
   function [ Q ] = lila_ormqrbz_v05_w02( m, n, k, i, j, mt, A, T, Q )
%
   Q(i:i+k-1,j:j+n-1) = zeros(k,n);
%
   V = tril( A( i:m,i:i+k-1 ), -1 ) + eye( m-i+1, k );
%
   H = ( eye( m-i+1,m-i+1 ) - V * ( T(i:i+k-1,i:i+k-1) * V' ) );
%
   Q( i:m,j:j+n-1 ) = H * Q( i:m,j:j+n-1 );
%
   end
