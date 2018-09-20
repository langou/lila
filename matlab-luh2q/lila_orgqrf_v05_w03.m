%
   function [ Q ] = lila_orgqrf_v05_w03( m, n, i, mt, A, T, Q )
%
   ilo = i;
   ihi = ilo + n - 1;
%
QQ = Q;
   QQ(ilo:m,ilo:ihi) = eye( size ( QQ(ilo:m,ilo:ihi) ) );
  for j = ihi:-1:ilo,
     [ QQ(j:m,j:ihi) ] = lapack_larfL( A(j:m,j), QQ(j:m,j:ihi) );
  end
%
%fprintf('cheating   %d\n',ihi:-1:ilo);


%
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
%fprintf('\n');
      while ( not_done == 1 )
%
%fprintf('broken   %d\n',jahi:-1:jalo);

        V = tril(A(ialo:iahi,jalo:jahi), -1) + eye(ml,vb);
%
        H = (eye(ml,ml) - V * ( T(itlo:ithi,jtlo:jthi) * V' ) );
%
        Q(ialo:iahi,ialo:ihi) = H*Q(ialo:iahi,ialo:ihi);
%       Q(ialo:iahi,jalo:jahi) = H*Q(ialo:iahi,jalo:jahi);
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
            itlo = mod(itlo-vb,mt);
            ithi = itlo+vb-1;
            jtlo = jtlo-vb;
%
         end
%
      end
%

   norm( QQ(ilo:m,ilo:ihi) - Q(ilo:m,ilo:ihi))
Q = QQ;


%fprintf('\n');
   end

