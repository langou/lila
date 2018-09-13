%
   function [ A, Q, TT, T ] = lila_geqrf_attempt_v05_w01_lvl1( m, nb, i, A, Q, TT, mt, T, nb_lvl2 );
%
   if ( m < i ) fprintf('m < i\n'); return; end
%
   nb_block = size(nb,2);
%
   ilo = zeros(1,nb_block);
   ilo(1) = 1;
   for ii = 2:nb_block,
      ilo(ii) = ilo(ii-1) + nb(ii-1);
   end
%
   ihi = zeros(1,nb_block);
   ihi(1) = nb(1);
   for ii = 2:nb_block,
      ihi(ii) = ihi(ii-1) + nb(ii);
   end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  this will go away at some point
   fprintf('\n')
   n = sum(nb_block);
   As = A;
%
   lda = -1;
   ldq = -1;
   ldt = -1;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%  GEQRF on the first block
%
   for j = 1:ihi(1), 
      [ A(j:m,j) ] = larfg( A(j:m,j) );
      [ A(j:m,j+1:ihi(1)) ] = larfL( A(j:m,j), A(j:m,j+1:ihi(1)) );
   end
%
   T(1:ihi(1),1:ihi(1)) = larft( A(1:m,1:ihi(1)) );
%
%  ORGQR on the first block
%
%  R = triu(A(1:ihi(1),1:ihi(1)));
   Q(1:m,ilo(1):ihi(1)) = eye( m, ihi(1) );
   for j = ihi(1):-1:1,
      [ Q(j:m,j:ihi(1)) ] = larfL( A(j:m,j), Q(j:m,j:ihi(1)) );
   end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
   if ( mt  > nb(1) ) vt = nb(1); else vt = mt; end;
%
   ilot = 1;
   ihit = vt;
   jlot = 1;
   jhit = vt;
%
   while ( nb(1) ~= jhit )
%  fprintf('jlot=%d, jhit=%d, vt=%d\n',jlot,jhit,vt);
   TT(ilot:ihit,jlot:jhit) = T(jlot:jhit,jlot:jhit);
   if ( jhit + mt > nb(1) ) vt = nb(1) - jhit; else vt = mt; end;
   ilot = 1;
   ihit = vt;
   jlot = jhit+1;
   jhit = jhit+vt;
   end
%  fprintf('jlot=%d, jhit=%d, vt=%d\n',jlot,jhit,vt);
   TT(ilot:ihit,jlot:jhit) = T(jlot:jhit,jlot:jhit);
   if ( jhit + mt > nb(1) ) vt = nb(1) - jhit; else vt = mt; end;
%
%  Second block
%
   [ A ] = lila_ormqrf_v05( m, nb(2), ihi(1), i, A, 1, 1, lda, A, 1, ilo(2), lda, TT, mt );
%
   [ A, Q, TT, T ] = lila_geqrf_attempt_v05_w01_lvl2( m, nb_lvl2{2}, ilo(2), A, Q, TT, mt, T );
%
   vb= mt - mod(ilo(2)-1,mt) ;
   if ( vb > nb(2) ), vb = nb(2); end
%
   itlo = mod(ilo(2)-1,mt)+1;
   ithi = itlo+vb-1;
   jtlo = ilo(2);
   jthi = jtlo+vb-1;
%
   TT(1:itlo-1,jtlo:jthi) = ...
       - ( TT(1:itlo-1,jtlo-(itlo)+1:jtlo-1) * A(jtlo:m,jtlo-(itlo)+1:jtlo-1)' )* ...
           ( (tril(A(jtlo:m,jtlo:jthi),-1) + eye(size((A(jtlo:m,jtlo:jthi)))) ) * TT(itlo:ithi,jtlo:jthi) );
%

%
   [ Q ] = lila_ormqrbz_v05( m, nb(2), ihi(1), i, A, 1, 1, lda, Q, 1, ilo(2), ldq, T, 1, 1, ldt );
%
