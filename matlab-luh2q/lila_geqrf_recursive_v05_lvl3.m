%
   function [ A, Q, TT ] = lila_geqrf_recursive_v05_lvl3( m, nb, i, A, Q, TT, mt )
%
   ilo = i;
   ihi = ilo + nb - 1;
%
   for j = ilo:ihi, 
      [ A(j:m,j) ] = larfg( A(j:m,j) );
      [ A(j:m,j+1:ihi) ] = larfL( A(j:m,j), A(j:m,j+1:ihi) );
   end
%
   vb= mt - mod(ilo-1,mt) ;
   if ( vb > nb ), vb = nb; end
%
   itlo = mod(ilo-1,mt)+1;
   ithi = itlo+vb-1;
   jtlo = ilo;
   jthi = jtlo+vb-1;
%
   TT(itlo:ithi,jtlo:jthi) = larft( A(jtlo:m,jtlo:jthi ) );
%
   itlo = mod(itlo+vb-1,mt)+1;
   jtlo = jtlo+vb;
   if( jtlo+mt-1 >= ihi) vb = ihi + 1 - jtlo; else, vb = mt; end
   ithi = itlo+vb-1;
   jthi = jtlo+vb-1;
%
   while (( itlo == 1 )&&(vb~=0))
%
      TT(itlo:ithi,jtlo:jthi) = larft( A(jtlo:m,jtlo:jthi ) );
%
      itlo = mod(itlo+vb-1,mt)+1;
      jtlo = jtlo+vb;
      if( jtlo+mt-1 >= ihi) vb = ihi + 1 - jtlo; else, vb = mt; end
      ithi = itlo+vb-1;
      jthi = jtlo+vb-1;
%
   end
%
   Q(ilo:m,ilo:ihi) = eye( size ( Q(ilo:m,ilo:ihi) ) );
   for j = ihi:-1:ilo,
      [ Q(j:m,j:ihi) ] = larfL( A(j:m,j), Q(j:m,j:ihi) );
   end
%
end
