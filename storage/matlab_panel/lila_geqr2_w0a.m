%
   function [ A, T, Q ] = lila_geqr2_w0a( m, n, i, A, T, Q )
%
      [ A ] = lapack_geqr2( m, n, i, A );   
%
      [ T ] = lapack_larft( A );   
%
      [ Q ] = lapack_orgqr2( m, n, i, A, Q );
%
   end
