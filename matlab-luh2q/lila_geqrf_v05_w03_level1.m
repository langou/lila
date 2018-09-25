%
   function [ A, T, Q ] = lila_geqrf_v05_w03_level1( n_lvl, i_lvl, nb_lvl, m, n, i, mt, A, T, Q )
%
   nb = nb_lvl(i_lvl);
   nb_block = ceil( n / nb );
%
   i_lvl = i_lvl + 1;
%
   lda = -1;
   ldq = -1;
   ldt = -1;
%
   if ( nb > n ) vb = n; else vb = nb; end;
%
   if (i_lvl == n_lvl+1)
   [ A, T, Q ] = lila_geqr2_v05_w03( m, vb, i, mt, A, T, Q );
   else
   [ A, T, Q ] = lila_geqrf_v05_w03_level1( n_lvl, i_lvl, nb_lvl, m, vb, i, mt, A, T, Q );
   end
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
   if (i_lvl == n_lvl+1)
   [ A, T, Q ] = lila_geqr2_v05_w03( m, vb, ilo, mt, A, T, Q );
   else
   [ A, T, Q ] = lila_geqrf_v05_w03_level1( n_lvl, i_lvl, nb_lvl, m, vb, ilo, mt, A, T, Q );
   end
%
   [ T ] = lila_larft_connect_v05_w03(  m, vb, ilo, mt, A, T  );   
%
   [ Q ] = lila_ormqrbz_v05_w03( m, vb, kb, i, ilo, mt, A, T, Q );
%
   end
%
