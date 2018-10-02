%
   function [ Q ] = lila_ormqrbz_v05_w02_flopsave( m, n, k, i, j, mt, A, T, Q )
%
%   Q(i:i+k-1,j:j+n-1) = zeros(k,n);
%   V = tril( A( i:m,i:i+k-1 ), -1 ) + eye( m-i+1, k );
%   H = ( eye( m-i+1,m-i+1 ) - V * ( T(i:i+k-1,i:i+k-1) * V' ) );
%   Q( i:m,j:j+n-1 ) = H * Q( i:m,j:j+n-1 );
%
%    Q(i:i+k-1,j:j+n-1) = zeros(k,n);
%    V = tril( A( i:m,i:i+k-1 ), -1 ) + eye( m-i+1, k );
%    work = V'  * Q( i:m,j:j+n-1 );
%    work = T(i:i+k-1,i:i+k-1) * work;
%    Q( i:m,j:j+n-1 ) = Q( i:m,j:j+n-1 ) - V * work;
%
%    Q(i:i+k-1,j:j+n-1) = zeros(k,n);
%    work = (tril( A( i:m,i:i+k-1 ), -1 ) + eye( m-i+1, k ))'  * Q( i:m,j:j+n-1 );
%    work = T(i:i+k-1,i:i+k-1) * work;
%    Q( i:m,j:j+n-1 ) = Q( i:m,j:j+n-1 ) - (tril( A( i:m,i:i+k-1 ), -1 ) + eye( m-i+1, k )) * work;
% 
%   work(1:k,1:n) = A( i+k:m,i:i+k-1 )' * Q( i+k:m,j:j+n-1 );
%   work(1:k,1:n) = T(i:i+k-1,i:i+k-1) * work(1:k,1:n);
%   Q( i:i+k-1,j:j+n-1 ) =  - (tril( A( i:i+k-1,i:i+k-1 ), -1 ) + eye( k, k )) * work(1:k,1:n);
%   Q( i+k:m,j:j+n-1 ) = Q( i+k:m,j:j+n-1 ) - A( i+k:m,i:i+k-1 ) * work(1:k,1:n);
%
%
%
    lda = -1;
    ldq = -1;
    ldt = -1;
    ldwork = -1;
%
    work(1:k,1:n) = randn(k,n); % (we put randn to simulate an allocation as opposed to putting zeros)
    [ work ] = blas_gemm ( 'T', 'N', k, n, m-k-i+1, ( +1.e+00 ), A, i+k, i, lda, Q, i+k, j, ldq, ( 0.e+00 ), work, 1, 1, ldwork );
    [ work ] = blas_trmm ( 'L', 'U', 'N', 'N', k, n, ( +1.e+00 ), T, i, i, ldt, work, 1, 1, ldwork );
    Q( i:i+k-1,j:j+n-1 ) = work(1:k,1:n);
    [ Q ] = blas_trmm ( 'L', 'L', 'N', 'U', k, n, ( -1.e+00 ), A, i, i, lda, Q, i, j, ldq );
    [ Q ] = blas_gemm ( 'N', 'N', m-k-i+1, n, k, ( -1.e+00 ), A, i+k, i, lda, work, 1, 1, ldwork, ( +1.e+00 ), Q, i+k, j, ldq );
%
   end
