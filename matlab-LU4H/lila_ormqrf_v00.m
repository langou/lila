   function [ B ] = lila_ormqrf_v0( m, n, k, A, ia, ja, lda, B, ib, jb, ldb )
%
      ialo = ia;
      iahi = ia+m-1;
      jalo = ja;
%
      iblo = ib;
      ibhi = ib+m-1;
      jblo = jb;
      jbhi = jb+n-1;
%
      for c = 1:k,
         [ B(iblo:ibhi,jblo:jbhi) ] = larfL( A(ialo:iahi,jalo), B(iblo:ibhi,jblo:jbhi) );
         iblo = iblo + 1;
         ialo = ialo + 1;
         jalo = jalo + 1;
      end
%
   end
