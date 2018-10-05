   function [ T ] = lila_larft_connect_v05_w02_flopsave( m, n, i, mt, A, T )
%
%   T(1:i-1,i:i+n-1) = ...
%       - ( T(1:i-1,1:i-1) * ( A(i:m,1:i-1)' * ...
%            (tril(A(i:m,i:i+n-1),-1) + eye(m-i+1,n) ) ) * T(i:i+n-1,i:i+n-1) );
% 
%   work(1:i-1,1:n) = A(i:m,1:i-1)' * ( tril(A(i:m,i:i+n-1),-1) + eye(m-i+1,n) );
%   T(1:i-1,i:i+n-1) = - ( T(1:i-1,1:i-1) * work(1:i-1,1:n) * T(i:i+n-1,i:i+n-1) );
% 
%   T(1:i-1,i:i+n-1) = A(i:i+n-1,1:i-1)';
%   T(1:i-1,i:i+n-1) = T(1:i-1,i:i+n-1) * ( tril(A(i:i+n-1,i:i+n-1),-1) + eye(n,n) );
%   T(1:i-1,i:i+n-1) = T(1:i-1,i:i+n-1) + A(i+n:m,1:i-1)' * A(i+n:m,i:i+n-1);
%   T(1:i-1,i:i+n-1) = - triu( T(1:i-1,1:i-1) ) * T(1:i-1,i:i+n-1);
%   T(1:i-1,i:i+n-1) = T(1:i-1,i:i+n-1) * triu( T(i:i+n-1,i:i+n-1) );
% 
    lda = -1;
    ldt = -1;
%
    T(1:i-1,i:i+n-1) = A(i:i+n-1,1:i-1)';
    [ T ] = blas_trmm ( 'R', 'L', 'N', 'U', i-1, n, (+1.0e+00), A, i, i, lda, T, 1, i, ldt );
    [ T ] = blas_gemm ( 'T', 'N', i-1, n, m-n-i+1, (+1.0e+00), A, i+n, 1, lda, A, i+n, i, lda, (+1.0e+00), T, 1, i, ldt );
    [ T ] = blas_trmm ( 'L', 'U', 'N', 'N', i-1, n, (-1.0e+00), T, 1, 1, ldt, T, 1, i, ldt );
    [ T ] = blas_trmm ( 'R', 'U', 'N', 'N', i-1, n, (+1.0e+00), T, i, i, ldt, T, 1, i, ldt );
% 
   end
