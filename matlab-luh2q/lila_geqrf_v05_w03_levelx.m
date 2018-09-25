%
   function [ A, T, Q ] = lila_geqrf_v05_w03_levelx( ii_lvl, nb_lvl, m, n, i, mt, A, T, Q )
%
   nb = nb_lvl(ii_lvl);
   nb_block = ceil( n / nb );
%
   lda = -1;
   ldq = -1;
   ldt = -1;
%
   if ( nb > n ) vb = n; else vb = nb; end;

   [ A, T, Q ] = lila_geqr2_v05_w03( m, vb, i, mt, A, T, Q );
%
   ilo = i;
   kb = 0;
%
   for j = 2:nb_block,
%
   ilo = ilo + vb;
   kb = kb + vb;
   if ( kb + nb > n ) vb = n - kb ; else vb = nb; end;
%
   [ A ] = lila_ormqrf_v05_w03( m, vb, kb, i, ilo, mt, A, T );
%
   [ A, T, Q ] = lila_geqr2_v05_w03( m, vb, ilo, mt, A, T, Q );
%
   [ T ] = lila_larft_connect_v05_w03(  m, vb, ilo, mt, A, T  );   
%
   [ Q ] = lila_ormqrbz_v05_w03( m, vb, kb, i, ilo, mt, A, T, Q );
%
   end
