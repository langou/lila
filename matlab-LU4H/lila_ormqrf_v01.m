   function [ B ] = lila_ormqrf_v01( m, n, k, A, ia, ja, lda, B, ib, jb, ldb, T )
%
      ialo = ia;
      iahi = ia+m-1;
      jalo = ja;
%
      iblo = ib;
      ibhi = ib+m-1;
      jblo = jb;
      jbhi = jb+k-1;
%
%      for c = 1:n,
%         [ B(iblo:ibhi,jblo:jbhi) ] = larfL( A(ialo:iahi,jalo), B(iblo:ibhi,jblo:jbhi) );
%         iblo = iblo + 1;
%         ialo = ialo + 1;
%         jalo = jalo + 1;
%      end
%
      V = tril(A(1:m,1:n), -1) + eye(m,n);
      H = (eye(m,m) - V * ( T(1:n,1:n)' * V' ) );
      B(1:m,jb:jb+k-1) = H*B(1:m,jb:jb+k-1);

   end
