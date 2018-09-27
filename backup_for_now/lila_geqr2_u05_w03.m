%
   function [ A, T ] = lila_geqr2_u05_w03( m, n, i, mt, A, T )
%
   [ A ] = lapack_geqr2( m, n, i, A );   
%
   [ T ] = lila_larft_u05_w03( m, n, i, mt, A, T );  
%
   end
