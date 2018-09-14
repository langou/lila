%
   function [ A, Q, TT ] = lila_geqrf_recursive_v05( m, n, i, A, Q, TT, mt )
%
if ( n < 10 )

   [ A, Q, TT ] = lila_geqrf_recursive_v05_lvl3( m, n, i, A, Q, TT, mt );

% if (n == 1)

   % [ Q(i:m,i) ] = A(i:m,i) / norm(  A(i:m,i) );
   % [ A(i:m,i) ] = larfg( A(i:m,i) );
   % if (A(i,i)<0), Q(i:m,i) = -Q(i:m,i); end;
   % it = mod(i,mt); if(it==0) it=mt; end;
   % TT(it,i) = larft( A(i:m,i) );

else
%
   nb(1) = ceil(n/2);
   nb(2) = n-nb(1);
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
   [ A, Q, TT ] = lila_geqrf_recursive_v05( m, nb(1), ilo(1), A, Q, TT, mt );
%
   [ A ] = lila_ormqrf_v05( m, nb(2), nb(1), ilo(1), A, 1, 1, lda, A, 1, ilo(2), lda, TT, mt );
%
   [ A, Q, TT ] = lila_geqrf_recursive_v05( m, nb(2), ilo(2), A, Q, TT, mt );
%
   [ TT ] = lila_larft_connect_v05( m, nb(2), ilo(2), A, TT, mt  );
%
   [ Q ] = lila_ormqrbz_v05( m, nb(2), nb(1), ilo(1), A, 1, 1, lda, Q, 1, ilo(2), ldq, TT, mt );
%
end
