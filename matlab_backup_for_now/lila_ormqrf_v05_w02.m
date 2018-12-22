   function [ A ] = lila_ormqrf_v05_w02( m, n, k, i, j, mt, A, T )
%
   V = tril(A(i:m,i:i+k-1), -1) + eye(m-i+1,k);
   H = (eye(m-i+1,m-i+1) - V * ( T(i:i+k-1,i:i+k-1)' * V' ) );
   A(i:m,j:j+n-1) = H*A(i:m,j:j+n-1);
%
   end
