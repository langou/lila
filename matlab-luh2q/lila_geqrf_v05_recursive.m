%
   function [ A, T, Q ] = lila_geqrf_v05_recursive( m, n, i, mt, A, T, Q )
%
if ( n < 10 )

   [ A, T, Q ] = lila_geqr2_v05( m, n, i, mt, A, T, Q );

% if (n == 1)

   % [ Q(i:m,i) ] = A(i:m,i) / norm(  A(i:m,i) );
   % [ A(i:m,i) ] = lapack_larfg( A(i:m,i) );
   % if (A(i,i)<0), Q(i:m,i) = -Q(i:m,i); end;
   % it = mod(i,mt); if(it==0) it=mt; end;
   % T(it,i) = lapack_larft( A(i:m,i) );

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
%   NOTE that I CHEAT:
%      I use larft(V) within some of the subroutines
%      This will need to change when we get the big T back
%
   [ A, T, Q ] = lila_geqrf_v05_recursive( m, nb(1), ilo(1), mt, A, T, Q );
%
%   [ A ] = lila_ormqrf_v05_w01( m, nb(2), nb(1), ilo(1), ilo(2), mt, A, T );
   [ A ] = lila_ormqrf_v05_w02( m, nb(2), nb(1), ilo(1), ilo(2), mt, A, T );
%   [ A ] = lila_ormqrf_v05_w03( m, nb(2), nb(1), ilo(1), ilo(2), mt, A, T );
%
   [ A, T, Q ] = lila_geqrf_v05_recursive( m, nb(2), ilo(2), mt, A, T, Q );
%
%   [ T ] = lila_larft_connect_v05_w01( m, nb(2), ilo(2), mt, A, T  );   %- DNW
   [ T ] = lila_larft_connect_v05_w02( m, nb(2), ilo(2), mt, A, T  );   
%   [ T ] = lila_larft_connect_v05_w03( m, nb(2), ilo(2), mt, A, T  );   %- DNW
%
%   [ Q ] = lila_ormqrbz_v05_w01( m, nb(2), nb(1), ilo(1), ilo(2), mt, A, T, Q );
   [ Q ] = lila_ormqrbz_v05_w02( m, nb(2), nb(1), ilo(1), ilo(2), mt, A, T, Q );
%   [ Q ] = lila_ormqrbz_v05_w03( m, nb(2), nb(1), ilo(1), ilo(2), mt, A, T, Q );   %- DNW
%
end
