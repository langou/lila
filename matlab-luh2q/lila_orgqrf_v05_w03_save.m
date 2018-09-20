%
   function [ Q ] = lila_orgqrf_v05_w03( m, n, i, mt, A, T, Q )
%
   ilo = i;
   ihi = ilo + n - 1;
%
   Q(ilo:m,ilo:ihi) = eye( size ( Q(ilo:m,ilo:ihi) ) );
   for j = ihi:-1:ilo,
      [ Q(j:m,j:ihi) ] = lapack_larfL( A(j:m,j), Q(j:m,j:ihi) );
   end
%
end



%
   function [ Q ] = lila_ormqrbz_v05_w03( m, n, k, i, j, mt, A, T, Q )
%
      vb = mod(i+k-1,mt); if ( vb == 0), vb = mt; end;
%
      if ( vb > k ), vb = k; end
%
      ialo = i+k-vb;
      iahi = m;
      jalo = i+k-vb;
      jahi = jalo+vb-1;
%
      iblo = i+k-vb;
      ibhi = m;
      jblo = j;
      jbhi = j+n-1;
%
      ithi = mod(i+k-1,mt); if (ithi == 0), ithi = mt; end
      itlo = ithi-vb+1;
      jtlo = i+k-vb;
      jthi = jtlo+vb-1;
%
      ml = m-i-k+1+vb ;
%
      not_done = 1;
%
      jj = 1;
%
      Q(i:i+k-1,jblo:jbhi) = zeros( size( Q(i:i+k-1,jblo:jbhi) ));
%
      while ( not_done == 1 )
%
        V = tril(A(ialo:iahi,jalo:jahi), -1) + eye(ml,vb);
%
        H = (eye(ml,ml) - V * ( T(itlo:ithi,jtlo:jthi) * V' ) );
%
        Q(iblo:ibhi,jblo:jbhi) = H*Q(iblo:ibhi,jblo:jbhi);
%
         if ( jj + vb - 1 == k ) 
%
            not_done = 0; 
%
         else if ( jj + vb - 1 > k ) 
%
             fprintf('defensive programming, we should never have been there, abort\n'); return;
%
         else
%
            jahi = jahi-vb;
            jthi = jthi-vb;
%
            jj = jj + vb;
%
            if ( ( jj+mt-1 ) <= k ), vb = mt; else vb = k-jj+1; end;
%
            ml = ml+vb;
%
            ialo = ialo-vb;
            iblo = iblo-vb;
            jalo = jalo-vb;
%
            itlo = mod(itlo-vb,mt);
            ithi = itlo+vb-1;
            jtlo = jtlo-vb;
%
         end
%
      end
%
   end

