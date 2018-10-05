
function [ B ] = blas_trmm ( side, uplo, transa, diag, m, n, alpha, A, ia, ja, lda, B, ib, jb, ldb )

   case_ = 'xxxx';
   if (( side == 'L' ) && ( uplo == 'L' ) && ( transa == 'T' ) && ( diag == 'U' )), case_ = 'LLTU'; end
   if (( side == 'L' ) && ( uplo == 'L' ) && ( transa == 'N' ) && ( diag == 'U' )), case_ = 'LLNU'; end
   if (( side == 'L' ) && ( uplo == 'U' ) && ( transa == 'T' ) && ( diag == 'N' )), case_ = 'LUTN'; end
   if (( side == 'L' ) && ( uplo == 'U' ) && ( transa == 'N' ) && ( diag == 'N' )), case_ = 'LUNN'; end
   if (( side == 'R' ) && ( uplo == 'U' ) && ( transa == 'N' ) && ( diag == 'N' )), case_ = 'RUNN'; end
   if (( side == 'R' ) && ( uplo == 'L' ) && ( transa == 'N' ) && ( diag == 'U' )), case_ = 'RLNU'; end

   switch (case_)
   case 'LLTU'
      B(ib:ib+m-1,jb:jb+n-1) = alpha * ((tril(A(ia:ia+m-1,ja:ja+m-1),-1)+eye(m,m))'*B(ib:ib+m-1,jb:jb+n-1));
   case 'LLNU'
      B(ib:ib+m-1,jb:jb+n-1) = alpha * ((tril(A(ia:ia+m-1,ja:ja+m-1),-1)+eye(m,m))*B(ib:ib+m-1,jb:jb+n-1));
   case 'LUTN'
      B(ib:ib+m-1,jb:jb+n-1) = alpha * (triu(A(ia:ia+m-1,ja:ja+m-1)))'*B(ib:ib+m-1,jb:jb+n-1);
   case 'LUNN'
      B(ib:ib+m-1,jb:jb+n-1) = alpha * (triu(A(ia:ia+m-1,ja:ja+m-1)))*B(ib:ib+m-1,jb:jb+n-1);
   case 'RUNN'
      B(ib:ib+m-1,jb:jb+n-1) = alpha * B(ib:ib+m-1,jb:jb+n-1) * (triu(A(ia:ia+n-1,ja:ja+n-1)));
   case 'RLNU'
      B(ib:ib+m-1,jb:jb+n-1) = alpha * B(ib:ib+m-1,jb:jb+n-1) * (tril(A(ia:ia+n-1,ja:ja+n-1),-1)+eye(n,n));
   otherwise
   fprintf('the side/uplo/transa/diag calling sequence of one blas_trmm is not correct\n');
   fprintf('or not yet implemented\n');
   end

end
