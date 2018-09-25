%
   function [ Q ] = lila_orgqrf_v05_w04( m, n, i, mt, A, T, Q )
%
   ilo = i;
   ihi = ilo + n - 1;
%
   Q(ilo:m,ilo:ihi) = eye( size ( Q(ilo:m,ilo:ihi) ) );
%
      vb = mod(i+n-1,mt); if ( vb == 0), vb = mt; end;
%
      if ( vb > n ), vb = n; end
%
      ialo = i+n-vb;
      iahi = m;
      jalo = i+n-vb;
      jahi = jalo+vb-1;
%
      ithi = mod(i+n-1,mt); if (ithi == 0), ithi = mt; end
      itlo = ithi-vb+1;
      jtlo = i+n-vb;
      jthi = jtlo+vb-1;
%
      ml = m-i-n+1+vb ;
%
      not_done = 1;
%
      jj = 1;
%
      while ( not_done == 1 )
%
        V = tril(A(ialo:iahi,jalo:jahi), -1) + eye(ml,vb);
%
        H = (eye(ml,ml) - V * ( T(itlo:ithi,jtlo:jthi) * V' ) );
%
        Q(ialo:m,ialo:ihi) = H*Q(ialo:m,ialo:ihi);
%
         if ( jj + vb - 1 == n ) 
%
            not_done = 0; 
%
         else if ( jj + vb - 1 > n ) 
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
            if ( ( jj+mt-1 ) <= n ), vb = mt; else vb = n-jj+1; end;
%
            ml = ml+vb;
%
            ialo = ialo-vb;
            jalo = jalo-vb;
%
            itlo = mod(itlo-vb,mt); if (itlo == 0), itlo = mt; end
            ithi = itlo+vb-1;
            jtlo = jtlo-vb;
%
         end
%
      end
%
   end

