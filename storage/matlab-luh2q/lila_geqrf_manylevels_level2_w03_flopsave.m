%
   function [ A, T, Q ] = lila_geqrf_manylevels_level2_w03_flopsave( m, n, i, mt, A, T, Q )
%
   global ii_lvl2;
   global nb_lvl2;
%
   nb = nb_lvl2{ii_lvl2};
%
   nb_block = size(nb,2);
%
   ilo = zeros(1,nb_block);
   ihi = zeros(1,nb_block);
   ilo(1) = i;
   ihi(1) = ilo(1) + nb(1) - 1;
   for ii = 2:nb_block,
      ilo(ii) = ilo(ii-1) + nb(ii-1);
      ihi(ii) = ihi(ii-1) + nb(ii);
   end
%
   lda = -1;
   ldq = -1;
   ldt = -1;
%
%  [ A, T, Q ] = lila_geqr2_v05_w03_a( m, nb(1), ilo(1), mt, A, T, Q );
   [ A, T, Q ] = lila_geqr2_v05_w03_b( m, nb(1), ilo(1), mt, A, T, Q );
%
   for k = 2:nb_block,
%
   [ A ] = lila_ormqrf_v05_w03_flopsave( m, nb(k), ihi(k-1)-ilo(1)+1, ilo(1), ilo(k), mt, A, T );
%
%  [ A, T, Q ] = lila_geqr2_v05_w03_a( m, nb(k), ilo(k), mt, A, T, Q );
   [ A, T, Q ] = lila_geqr2_v05_w03_b( m, nb(k), ilo(k), mt, A, T, Q );
%
   [ T ] = lila_larft_connect_v05_w03_flopsave(  m, nb(k), ilo(k), mt, A, T  );   
%
   [ Q ] = lila_ormqrbz_v05_w03_flopsave( m, nb(k), ihi(k-1)-ilo(1)+1, ilo(1), ilo(k), mt, A, T, Q );
%
   end
