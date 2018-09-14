%
   function [ A, Q, TT ] = lila_geqrf_attempt_v05_w01_lvl2( m, n, i, A, Q, TT, mt )
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  this will go away at some point
   fprintf('\n')
   n = sum(nb_block);
   if ( m < i+n ) fprintf('m < i+n\n'); return; end
%
   lda = -1;
   ldq = -1;
   ldt = -1;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
   [ A, Q, TT ] = lila_geqrf_attempt_v05_w01_lvl3( m, nb(1), ilo(1), A, Q, TT, mt );
%
   for k = 2:nb_block,
%
%%%%%%%%%%%%%%%
%
   [ A ] = lila_ormqrf_v05( m, nb(k), ihi(k-1)-ilo(1)+1, ilo(1), A, 1, 1, lda, A, 1, ilo(k), lda, TT, mt );
%
   [ A, Q, TT ] = lila_geqrf_attempt_v05_w01_lvl3( m, nb(k), ilo(k), A, Q, TT, mt );
%
   [ TT ] = lila_larft_connect_v05( m, nb(k), ilo(k), A, TT, mt  );
%
   [ Q ] = lila_ormqrbz_v05( m, nb(k), ihi(k-1)-ilo(1)+1, ilo(1), A, 1, 1, lda, Q, 1, ilo(k), ldq, TT, mt );
%
   end
