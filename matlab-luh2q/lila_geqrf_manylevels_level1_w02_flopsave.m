%
   function [ A, T, Q ] = lila_geqrf_manylevels_level1_w02_flopsave( m, n, i, mt, A, T, Q )
%
   global nb_lvl1;
   global ii_lvl2;
   global nb_lvl2;
%
   nb = nb_lvl1;
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
   ii_lvl2 = 1;
   [ A, T, Q ] = lila_geqrf_manylevels_level2_w02_flopsave( m, nb(1), ilo(1), mt, A, T, Q );
%
   for k = 2:nb_block,
%
   ii_lvl2 = k;
%
   [ A ] = lila_ormqrf_v05_w02_flopsave( m, nb(k), ihi(k-1)-ilo(1)+1, ilo(1), ilo(k), mt, A, T );
%
   [ A, T, Q ] = lila_geqrf_manylevels_level2_w02_flopsave( m, nb(k), ilo(k), mt, A, T, Q );
%
   [ T ] = lila_larft_connect_v05_w02( m, nb(k), ilo(k), mt, A, T );  
%
   [ Q ] = lila_ormqrbz_v05_w02( m, nb(k), ihi(k-1)-ilo(1)+1, ilo(1), ilo(k), mt, A, T, Q );
%
   end
