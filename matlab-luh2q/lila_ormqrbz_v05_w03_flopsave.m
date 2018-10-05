%
   function [ Q ] = lila_ormqrbz_v05_w03_flopsave( m, n, k, i, j, mt, A, T, Q )
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
%        V = tril(A(ialo:iahi,jalo:jahi), -1) + eye(ml,vb);
%        H = (eye(ml,ml) - V * ( T(itlo:ithi,jtlo:jthi) * V' ) );
%        Q(iblo:ibhi,jblo:jbhi) = H*Q(iblo:ibhi,jblo:jbhi);
%
%        V = tril(A(ialo:iahi,jalo:jahi), -1) + eye(ml,vb);
%        Q(iblo:ibhi,jblo:jbhi) = (Q(iblo:ibhi,jblo:jbhi) - V * ( T(itlo:ithi,jtlo:jthi) * V' ) * Q(iblo:ibhi,jblo:jbhi) );
%
%        V = tril(A(ialo:iahi,jalo:jahi), -1) + eye(ml,vb);
%        work = V' * Q(iblo:ibhi,jblo:jbhi);
%        work = T(itlo:ithi,jtlo:jthi) * work;
%        work = V * work;
%        Q(iblo:ibhi,jblo:jbhi) = (Q(iblo:ibhi,jblo:jbhi) - work );
%
%        work = zeros(ml,n);
%        work(1:vb,1:n) = (tril(A(ialo:iahi,jalo:jahi), -1) + eye(ml,vb))' * Q(iblo:ibhi,jblo:jbhi);
%        work(1:vb,1:n) = T(itlo:ithi,jtlo:jthi) * work(1:vb,1:n);
%        work(1:ml,1:n) = (tril(A(ialo:iahi,jalo:jahi), -1) + eye(ml,vb)) * work(1:vb,1:n);
%        Q(iblo:ibhi,jblo:jbhi) = Q(iblo:ibhi,jblo:jbhi) - work(1:ml,1:n);
%%%%%%
% for reference
%%%%%%
%        work = zeros(ml,n);
%        work(1:vb,1:n) = (tril(A(ialo:jahi,jalo:jahi), -1) + eye(vb,vb))' * Q(iblo:iblo+vb-1,jblo:jbhi);
%        work(1:vb,1:n) = A(jahi+1:iahi,jalo:jahi)' * Q(iblo+vb:ibhi,jblo:jbhi);
%        work(1:vb,1:n) = T(itlo:ithi,jtlo:jthi) * work(1:vb,1:n);
%%%%%%%
%  Stuck trying to break this one apart
%        work(1:ml,1:n) = (tril(A(ialo:iahi,jalo:jahi), -1) + eye(ml,vb)) * work(1:vb,1:n);
%        work(1:vb,1:vb) = (tril(A(ialo:ialo+vb-1,jalo:jahi), -1) + eye(vb,vb)) * work(1:vb,1:vb);
%        work(vb+1:ml,vb+1:n) = A(ialo+vb:iahi,jalo:jahi) * work(1:vb,vb+1:n);
%%%%%%%
%        Q(iblo:ibhi,jblo:jbhi) = Q(iblo:ibhi,jblo:jbhi) - work(1:ml,1:n);
%%%%%
%
        lda = -1;
        ldb = -1;
        ldwork = -1;
%
        work = zeros(ml,n);
        [ work ] = blas_trmm( 'L', 'L', 'T', 'U', vb, n, (+1.0e00), A, ialo, jalo, lda, Q, iblo, jblo, ldb );
        work(1:vb,1:n) = A(jahi+1:iahi,jalo:jahi)' * Q(iblo+vb:ibhi,jblo:jbhi);
        work(1:vb,1:n) = T(itlo:ithi,jtlo:jthi) * work(1:vb,1:n);
        work(1:ml,1:n) = (tril(A(ialo:iahi,jalo:jahi), -1) + eye(ml,vb)) * work(1:vb,1:n);
        Q(iblo:ibhi,jblo:jbhi) = Q(iblo:ibhi,jblo:jbhi) - work(1:ml,1:n);
%
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
            itlo = mod(itlo-vb,mt); if (itlo == 0), itlo = mt; end
            ithi = itlo+vb-1;
            jtlo = jtlo-vb;
%
         end
%
      end
%
   end
