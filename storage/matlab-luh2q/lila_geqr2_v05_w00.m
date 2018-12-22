%
   function [ A, T, Q ] = lila_geqr2_v05_w00( m, n, i, mt, A, T, Q )
%
   [ A ] = lapack_geqr2( m, n, i, A );   
%
   [ T ] = rand( 1, n );
%
   [ Q ] = lila_orgqrf_v05_w00( m, n, i, mt, A, T, Q );
%
end
