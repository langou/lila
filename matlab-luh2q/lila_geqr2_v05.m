%
   function [ A, T, Q ] = lila_geqr2_v05( m, n, i, mt, A, T, Q )
%
   ilo = i;
   ihi = ilo + n - 1;
%
   [ A ] = lapack_geqr2( m, n, i, A );
%
   [ T ] = lila_larft_v05( m, n, i, mt, A, T );
%
   [ Q ] = lapack_orgqr( m, n, i, A, Q );
%
end
