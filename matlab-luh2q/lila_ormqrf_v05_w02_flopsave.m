   function [ A ] = lila_ormqrf_v05_w02_flopsave( m, n, k, i, j, mt, A, T )
%
%  V = tril(A(i:m,i:i+k-1), -1) + eye(m-i+1,k);
%  H = (eye(m-i+1,m-i+1) - V * ( T(i:i+k-1,i:i+k-1)' * V' ) );
%  A(i:m,j:j+n-1) = H*A(i:m,j:j+n-1);
%
%  V = tril(A(i:m,i:i+k-1), -1) + eye(m-i+1,k);
%  A(i:m,j:j+n-1) = A(i:m,j:j+n-1) - V * ( T(i:i+k-1,i:i+k-1)' * ( V'  * A(i:m,j:j+n-1) ) );
%
%  V = tril(A(i:m,i:i+k-1), -1) + eye(m-i+1,k);
%  work = V' * A(i:m,j:j+n-1);
%  work = T(i:i+k-1,i:i+k-1)' * work;
%  A(i:m,j:j+n-1) = A(i:m,j:j+n-1) - V * work;
%
%  work = ( tril(A(i:m,i:i+k-1), -1) + eye(m-i+1,k) )' * A(i:m,j:j+n-1);
%  work = T(i:i+k-1,i:i+k-1)' * work;
%  A(i:m,j:j+n-1) = A(i:m,j:j+n-1) - ( tril(A(i:m,i:i+k-1), -1) + eye(m-i+1,k) ) * work;
%
%  work(1:k,1:n) = ( tril(A(i:i+k-1,i:i+k-1), -1) + eye(k,k) )' * A(i:i+k-1,j:j+n-1);
%  work(1:k,1:n) = work(1:k,1:n) + ( A(i+k:m,i:i+k-1) )' * A(i+k:m,j:j+n-1);
%  work(1:k,1:n) = T(i:i+k-1,i:i+k-1)' * work(1:k,1:n);
%  A(i:i+k-1,j:j+n-1) = A(i:i+k-1,j:j+n-1) - ( tril(A(i:i+k-1,i:i+k-1), -1) + eye(k,k) ) * work(1:k,1:n);
%  A(i+k:m,j:j+n-1) = A(i+k:m,j:j+n-1) - ( A(i+k:m,i:i+k-1) ) * work(1:k,1:n);
%
%
   lda = -1; 
   ldb = -1;
   ldt = -1;
   ldwork = -1;
%  allocate work of size k-by-n
%
   work = zeros(k,n);
   work(1:k,1:n) = A(i:i+k-1,j:j+n-1);
   [ work ] = blas_trmm ( 'L', 'L', 'T', 'U', k, n, (+1.0e+00), A, i, j, lda, work, 1, 1, ldb );
   [ work ] = blas_gemm ( 'T', 'N', k, n, m-k-i+1, ( +1.e+00 ), A, i+k, i, lda, A, i+k, j, lda, ( +1.e+00 ), work, 1, 1, ldwork );
   [ work ] = blas_trmm ( 'L', 'U', 'T', 'N', k, n, (+1.0e+00), T, i, i, ldt, work, 1, 1, ldwork );
   [ A ] = blas_gemm ( 'N', 'N', m-k-i+1, n, k, ( -1.e+00 ), A, i+k, i, lda, work, 1, 1, ldwork, ( +1.e+00 ), A, i+k, j, ldwork );
   [ work ] = blas_trmm ( 'L', 'L', 'N', 'U', k, n, (+1.0e+00), A, i, i, ldt, work, 1, 1, ldwork );
   A(i:i+k-1,j:j+n-1) = A(i:i+k-1,j:j+n-1) - work(1:k,1:n);
%
   end
