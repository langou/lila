%
   function [ A, T ] = lila_geqrf_u05_w03_level1( m, n, i, mt, A, T )
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
   [ A, T ] = lila_geqrf_u05_w03_level2( m, nb(1), ilo(1), mt, A, T );
%
   for k = 2:nb_block,
%
   ii_lvl2 = k;
%
   [ A ] = lila_ormqrf_u05_w03( m, nb(k), ihi(k-1)-ilo(1)+1, ilo(1), ilo(k), mt, A, T );
%
   [ A, T ] = lila_geqrf_u05_w03_level2( m, nb(k), ilo(k), mt, A, T );
%
   [ T ] = lila_larft_connect_u05_w03( m, nb(k), ilo(k), mt, A, T  );  
%
   end
