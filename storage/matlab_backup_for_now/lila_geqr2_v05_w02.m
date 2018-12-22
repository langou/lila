%
   function [ A, T, Q ] = lila_geqr2_v05_w02( m, n, i, mt, A, T, Q )
%
   [ A ] = lapack_geqr2( m, n, i, A );   
%
   [ T ] = lila_larft_v05_w02( m, n, i, mt, A, T );   
%
   [ Q ] = lila_orgqrf_v05_w02( m, n, i, mt, A, T, Q );
%
end
