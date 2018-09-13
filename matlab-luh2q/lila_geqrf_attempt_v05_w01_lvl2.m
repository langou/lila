%
   function [ A, Q, TT, T ] = lila_geqrf_attempt_v05_w01_lvl2( m, nb, i, A, Q, TT, mt, T );
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
%  GEQRF on the first block
%
   for j = ilo(1):ihi(1), 
      [ A(j:m,j) ] = larfg( A(j:m,j) );
      [ A(j:m,j+1:ihi(1)) ] = larfL( A(j:m,j), A(j:m,j+1:ihi(1)) );
   end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
   vb= mt - mod(ilo(1)-1,mt) ;
   if ( vb > nb(1) ), vb = nb(1); end
%
   itlo = mod(ilo(1)-1,mt)+1;
   ithi = itlo+vb-1;
   jtlo = ilo(1);
   jthi = jtlo+vb-1;
%
%  fprintf('itlo=%d, ithi=%d, jtlo=%d, jthi=%d, vb=%d\n',itlo,ithi,jtlo,jthi,vb);
   TT(itlo:ithi,jtlo:jthi) = larft( A(jtlo:m,jtlo:jthi ) );
%
   itlo = mod(itlo+vb-1,mt)+1;
   jtlo = jtlo+vb;
   if( jtlo+mt-1 >= ihi(1)) vb = ihi(1) + 1 - jtlo; else, vb = mt; end
   ithi = itlo+vb-1;
   jthi = jtlo+vb-1;
%
   while (( itlo == 1 )&&(vb~=0))
%
%     fprintf('itlo=%d, ithi=%d, jtlo=%d, jthi=%d, vb=%d\n',itlo,ithi,jtlo,jthi,vb);
      TT(itlo:ithi,jtlo:jthi) = larft( A(jtlo:m,jtlo:jthi ) );
%
      itlo = mod(itlo+vb-1,mt)+1;
      jtlo = jtlo+vb;
      if( jtlo+mt-1 >= ihi(1)) vb = ihi(1) + 1 - jtlo; else, vb = mt; end
      ithi = itlo+vb-1;
      jthi = jtlo+vb-1;
%
   end

%
   Q(ilo(1):m,ilo(1):ihi(1)) = eye( size ( Q(ilo(1):m,ilo(1):ihi(1)) ) );
   for j = ihi(1):-1:ilo(1),
      [ Q(j:m,j:ihi(1)) ] = larfL( A(j:m,j), Q(j:m,j:ihi(1)) );
   end

%%%%%%%%%%%%%%%
%
   [ A ] = lila_ormqrf_v05( m, nb(2), ihi(1)-ilo(1)+1, ilo(1), A, 1, 1, lda, A, 1, ilo(2), lda, TT, mt );
%
   for j = ilo(2):ihi(2), 
      [ A(j:m,j) ] = larfg( A(j:m,j) );
      [ A(j:m,j+1:ihi(2)) ] = larfL( A(j:m,j), A(j:m,j+1:ihi(2)) );
   end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
   vb= mt - mod(ilo(2)-1,mt) ;
   if ( vb > nb(2) ), vb = nb(2); end
%
   itlo = mod(ilo(2)-1,mt)+1;
   ithi = itlo+vb-1;
   jtlo = ilo(2);
   jthi = jtlo+vb-1;
%
%  fprintf('itlo=%d, ithi=%d, jtlo=%d, jthi=%d, vb=%d\n',itlo,ithi,jtlo,jthi,vb);
   TT(itlo:ithi,jtlo:jthi) = larft( A(jtlo:m,jtlo:jthi ) );
%
   itlo = mod(itlo+vb-1,mt)+1;
   jtlo = jtlo+vb;
   if( jtlo+mt-1 >= ihi(2)) vb = ihi(2) + 1 - jtlo; else, vb = mt; end
   ithi = itlo+vb-1;
   jthi = jtlo+vb-1;
%
   while (( itlo == 1 )&&(vb~=0))
%
%     fprintf('itlo=%d, ithi=%d, jtlo=%d, jthi=%d, vb=%d\n',itlo,ithi,jtlo,jthi,vb);
      TT(itlo:ithi,jtlo:jthi) = larft( A(jtlo:m,jtlo:jthi ) );
%
      itlo = mod(itlo+vb-1,mt)+1;
      jtlo = jtlo+vb;
      if( jtlo+mt-1 >= ihi(2)) vb = ihi(2) + 1 - jtlo; else, vb = mt; end
      ithi = itlo+vb-1;
      jthi = jtlo+vb-1;
%
   end





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
   Q(ilo(2):m,ilo(2):ihi(2)) = A(ilo(2):m,ilo(2):ihi(2));
   [ Q ] = lila_orgqr_v05( m, nb(2), ilo(2), Q, 1, 1, ldq, T, 1, 1, ldt );
%
   [ Q ] = lila_ormqrbz_v05( m, nb(2), ihi(1)-ilo(1)+1, ilo(1), A, 1, 1, lda, Q, 1, ilo(2), ldq, T, 1, 1, ldt );
%
