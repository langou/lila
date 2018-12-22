   function [ B ] = lila_ormqrf_v01( m, n, k, A, ia, ja, lda, B, ib, jb, ldb, T, it, jt, ldt )
%
      ialo = ia;
      iahi = ia+m-1;
      jalo = ja;
      jahi = ja+k-1;
%
      iblo = ib;
      ibhi = ib+m-1;
      jblo = jb;
      jbhi = jb+n-1;
%
      itlo = it;
      ithi = it+k-1;
      jtlo = jt;
      jthi = jt+k-1;
%
      V = tril(A(ialo:iahi,jalo:jahi), -1) + eye(m,k);
      H = (eye(m,m) - V * ( T(itlo:ithi,jtlo:jthi)' * V' ) );
      B(iblo:ibhi,jblo:jbhi) = H*B(iblo:ibhi,jblo:jbhi);
%
   end
