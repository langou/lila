   function [ TT ] = lila_larft_v05( m, n, i, mt, A, TT )
%
   ilo = i;
   ihi = ilo + n - 1;
%
   vb= mt - mod(ilo-1,mt) ;
   if ( vb > n ), vb = n; end
%
   itlo = mod(ilo-1,mt)+1;
   ithi = itlo+vb-1;
   jtlo = ilo;
   jthi = jtlo+vb-1;
%
   TT(itlo:ithi,jtlo:jthi) = lapack_larft( A(jtlo:m,jtlo:jthi ) );
%
   itlo = mod(itlo+vb-1,mt)+1;
   jtlo = jtlo+vb;
   if( jtlo+mt-1 >= ihi) vb = ihi + 1 - jtlo; else, vb = mt; end
   ithi = itlo+vb-1;
   jthi = jtlo+vb-1;
%
   while (( itlo == 1 )&&(vb~=0))
%
      TT(itlo:ithi,jtlo:jthi) = lapack_larft( A(jtlo:m,jtlo:jthi ) );
%
      itlo = mod(itlo+vb-1,mt)+1;
      jtlo = jtlo+vb;
      if( jtlo+mt-1 >= ihi) vb = ihi + 1 - jtlo; else, vb = mt; end
      ithi = itlo+vb-1;
      jthi = jtlo+vb-1;
%
   end
