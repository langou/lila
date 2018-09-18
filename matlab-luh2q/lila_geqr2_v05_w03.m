%
   function [ A, T, Q ] = lila_geqr2_v05_w03( m, n, i, mt, A, T, Q )
%
   ilo = i;
   ihi = ilo + n - 1;
%
   [ A ] = lapack_geqr2( m, n, i, A );   
%
%   [ T ] = lila_larft_v05_w01( m, n, i, mt, A, T );   
   [ T ] = lila_larft_v05_w02( m, n, i, mt, A, T );   
%   [ T ] = lila_larft_v05_w03( m, n, i, mt, A, T );   
%
   [ Q ] = lila_orgqrf_v05_w01( m, n, i, mt, A, T, Q );
%   [ Q ] = lila_orgqrf_v05_w02( m, n, i, mt, A, T, Q );   %- DNW 
%   [ Q ] = lila_orgqrf_v05_w03( m, n, i, mt, A, T, Q );  %- DNW
%
end
