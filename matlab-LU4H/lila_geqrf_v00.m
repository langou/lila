   function [ A ] = lila_geqrf_v0( m, n, A, ia, ja, lda )
%
      ilo = ia;
      ihi = ia+m-1;
      jlo = ja;
      jhi = ja+n-1;
%
      for c = 1:n,
         [ A( ilo:ihi, jlo ) ] = larfg( A( ilo:ihi, jlo ) );
         [ A( ilo:ihi, jlo+1:jhi ) ] = larfL( A( ilo:ihi, jlo ), A( ilo:ihi, jlo+1:jhi ) );
         ilo = ilo + 1;
         jlo = jlo + 1;
      end
%
   end
