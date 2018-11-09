%
   function [ A, T, Q, D ] = lila_orghr_w0b_panel( m, n, i, j, mt, A, T, Q, D )
%
      [ A, T ] = lila_ormhr_w0b( m, n, i, j, mt, A, T, Q, D );
%
      [ A, T, Q, D ] = lila_orgh2_w0b_panel( m, n, j, mt, A, T, Q, D );
%
      [ T ] = lila_larft_w0b_connect( n, i, j, mt, A, T );
%
   end
