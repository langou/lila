   function [ Q ] = lila_orgqrf_v05_w02( m, n, i, mt, A, T, Q )
%
   Q(i:m,i:i+n-1) = eye( m-i+1, n );
   V = tril(A(i:m,i:i+n-1), -1) + eye(m-i+1, n);
   H = (eye( m-i+1, m-i+1 ) - V * ( T(i:i+n-1,i:i+n-1) * V' ) );
   Q(i:m,i:i+n-1) = H*Q(i:m,i:i+n-1);
%
end
