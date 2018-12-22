%
   function [ Q ] = lila_orgqrf_v05_w01( m, n, i, mt, A, T, Q )
%
   Q(i:m,i:i+n-1) = eye( size ( Q(i:m,i:i+n-1) ) );
   for j = i+n-1:-1:i,
%
         V = tril(A(j:m,j), -1) + eye(size(A(j:m,1)));
%
         H = (eye(m-j+1,m-j+1) - V * ( T(1,j) * V' ) );
%
         Q(j:m,j:i+n-1) = H * Q(j:m,j:i+n-1);
%
   end
