   function [ A ] = lila_ormqrf_v05_w01_flopsave( m, n, k, i, j, mt, A, T )
%
%  for ii = i:i+k-1, 
%  V = tril(A(ii:m,ii), -1) + eye(m-ii+1,1);
%  H = (eye(m-ii+1,m-ii+1) - V * ( T(1,ii) * V' ) );
%  A(ii:m,j:j+n-1) = H*A(ii:m,j:j+n-1);
%  end
%
%  for ii = i:i+k-1, 
%  V = tril(A(ii:m,ii), -1) + eye(m-ii+1,1);
%  A(ii:m,j:j+n-1) = A(ii:m,j:j+n-1) - V *  T(1,ii) * ( V' * A(ii:m,j:j+n-1) );
%  end
%
%  for ii = i:i+k-1, 
%  V = tril(A(ii:m,ii), -1) + eye(m-ii+1,1);
%  work(1,1:n) = V' * A(ii:m,j:j+n-1);
%  work(1,1:n) = T(1,ii) * work(1,1:n);
%  A(ii:m,j:j+n-1) = A(ii:m,j:j+n-1) - V * work(1,1:n);
%  end
%
%   for ii = i:i+k-1, 
%   work(1,1:n) = (tril(A(ii:m,ii), -1) + eye(m-ii+1,1))' * A(ii:m,j:j+n-1);
%   work(1,1:n) = T(1,ii) * work(1,1:n);
%   A(ii:m,j:j+n-1) = A(ii:m,j:j+n-1) - ( tril(A(ii:m,ii), -1) + eye(m-ii+1,1) ) * work(1,1:n);
%   end
%
%   for ii = i:i+k-1, 
%   work(1,1:n) = [ 1; A(ii+1:m,ii) ]' * A(ii:m,j:j+n-1);
%   work(1,1:n) = T(1,ii) * work(1,1:n);
%   A(ii:m,j:j+n-1) = A(ii:m,j:j+n-1) - [ 1; A(ii+1:m,ii) ] * work(1,1:n);
%   end
%
%   for ii = i:i+k-1, 
%   work(1,1:n) = A(ii,j:j+n-1);
%   work(1,1:n) = work(1,1:n) + [ A(ii+1:m,ii) ]' * A(ii+1:m,j:j+n-1);
%   work(1,1:n) = T(1,ii) * work(1,1:n);
%   A(ii,j:j+n-1) = A(ii,j:j+n-1) - work(1,1:n);
%   A(ii+1:m,j:j+n-1) = A(ii+1:m,j:j+n-1) - [ A(ii+1:m,ii) ] * work(1,1:n);
%   end
%
   lda = -1;
   ldwork = -1;
    for ii = i:i+k-1, 
    work(1,1:n) = A(ii,j:j+n-1);
    work(1,1:n) = ( work(1,1:n)' + A(ii+1:m,j:j+n-1)' * ( A(ii+1:m,ii) ) )';
%    [ work ] = blas_gemm( 'T', 'N', 1, n, m, (+1.0e+00), A, ii+1, j, lda, A, ii+1, ii, lda, (+1.0e+00), work, 1, 1, ldwork )
    work(1,1:n) = T(1,ii) * work(1,1:n);
    A(ii,j:j+n-1) = A(ii,j:j+n-1) - work(1,1:n);
%    A(ii+1:m,j:j+n-1) = A(ii+1:m,j:j+n-1) - [ A(ii+1:m,ii) ] * work(1,1:n);
    [ A ] = blas_gemm( 'N', 'N', m-ii, n, 1, (-1.0e+00), A, ii+1, ii, lda, work, 1, 1, ldwork, (+1.0e+00), A, ii+1, j, lda );

    end
%
