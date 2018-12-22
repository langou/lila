   function [ B ] = lila_ormqrbz_v00( m, n, k, A, ia, ja, lda, B, ib, jb, ldb )
%
      ialo = ia+k-1;
      iahi = ia+m-1;
      jalo = ja+k-1;
%
      iblo = ib+k-1;
      ibhi = ib+m-1;
      jblo = jb;
      jbhi = jb+n-1;
%
      B(iblo-k+1:ibhi,jblo:jbhi) = [ zeros(k-1,n) ; B(iblo:ibhi,jblo:jbhi)];
%
      for c = k:-1:1,
         [ B(iblo:ibhi,jblo:jbhi) ] = larfL( A(ialo:iahi,jalo), B(iblo:ibhi,jblo:jbhi) );
         iblo = iblo - 1;
         ialo = ialo - 1;
         jalo = jalo - 1;
      end
%
   end
