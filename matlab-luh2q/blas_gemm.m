
function [ C ] = blas_gemm ( transa, transb, m, n, k, alpha, A, ia, ja, lda, B, ib, jb, ldb, beta, C, ic, jc, ldc )

   case_ = 'xx';
   if (( transa == 'T' ) && ( transb == 'N' )), case_ = 'TN'; end
   if (( transa == 'N' ) && ( transb == 'N' )), case_ = 'NN'; end

   switch (case_)
   case 'NN'
      C(ic:ic+m-1,jc:jc+n-1) = beta * C(ic:ic+m-1,jc:jc+n-1) + alpha * A(ia:ia+m-1,ja:ja+k-1) * B(ib:ib+k-1,jb:jb+n-1);
   case 'TN'
      C(ic:ic+m-1,jc:jc+n-1) = beta * C(ic:ic+m-1,jc:jc+n-1) + alpha * A(ia:ia+k-1,ja:ja+m-1)' * B(ib:ib+k-1,jb:jb+n-1);
   otherwise
   fprintf('the side/uplo/transa/diag calling sequence of one blas_trmm is not correct\n');
   fprintf('or not yet implemented\n');
   end

end
