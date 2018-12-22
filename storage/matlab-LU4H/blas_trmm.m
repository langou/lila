
function [ B ] = blas_trmm ( side, uplo, transa, diag, m, n, alpha, A, ia, ja, lda, B, ib, jb, ldb )

   case_ = 'xxxx';
   if (( side == 'L' ) && ( uplo == 'L' ) && ( transa == 'T' ) && ( diag == 'U' )), case_ = 'LLTU'; end
   if (( side == 'L' ) && ( uplo == 'L' ) && ( transa == 'N' ) && ( diag == 'U' )), case_ = 'LLNU'; end
   if (( side == 'L' ) && ( uplo == 'U' ) && ( transa == 'T' ) && ( diag == 'N' )), case_ = 'LUTN'; end

   switch (case_)
   case 'LLTU'
      B(ib:ib+m-1,jb:jb+n-1) = alpha * ((tril(A(ia:ia+m-1,ia:ia+m-1),-1)+eye(m,m))'*B(ib:ib+m-1,jb:jb+n-1));
   case 'LLNU'
      B(ib:ib+m-1,jb:jb+n-1) = alpha * ((tril(A(ia:ia+m-1,ia:ia+m-1),-1)+eye(m,m))*B(ib:ib+m-1,jb:jb+n-1));
   case 'LUTN'
      B(ib:ib+m-1,jb:jb+n-1) = alpha * (triu(A(ia:ia+m-1,ia:ia+m-1)))'*B(ib:ib+m-1,jb:jb+n-1);
   otherwise
   fprintf('the side/uplo/transa/diag calling sequence of one blas_trmm is not correct\n');
   fprintf('or not yet implemented\n');
   end


end
