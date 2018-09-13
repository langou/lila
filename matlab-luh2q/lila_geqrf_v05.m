   function [ A ] = lila_geqrf_v05( m, n, j, A, ia, ja, lda, T, mt )
%
      ilo = ia+j-1;
      ihi = m;
      jlo = ja+j-1;
      jhi = ja+j+n-1-1;
%
      for c = 1:n,
         [ A( ilo:ihi, jlo ) ] = larfg( A( ilo:ihi, jlo ) );
         [ A( ilo:ihi, jlo+1:jhi ) ] = larfL( A( ilo:ihi, jlo ), A( ilo:ihi, jlo+1:jhi ) );
         ilo = ilo + 1;
         jlo = jlo + 1;
      end
%
   end
