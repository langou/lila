   function [ B ] = lila_ormqrbz_v0( m, n, k, A, ia, ja, lda, B, ib, jb, ldb )
%
      ialo = ia+n-1;
      iahi = ia+m-1;
      jalo = ja+n-1;
%
      iblo = ib+n-1;
      ibhi = ib+m-1;
      jblo = jb;
      jbhi = jb+k-1;
%
      B(iblo-n+1:ibhi,jblo:jbhi) = [ zeros(n-1,k) ; B(iblo:ibhi,jblo:jbhi)];
%
      for c = n:-1:1,
         [ B(iblo:ibhi,jblo:jbhi) ] = larfL( A(ialo:iahi,jalo), B(iblo:ibhi,jblo:jbhi) );
         iblo = iblo - 1;
         ialo = ialo - 1;
         jalo = jalo - 1;
      end
%
   end
