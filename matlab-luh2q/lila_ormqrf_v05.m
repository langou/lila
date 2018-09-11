   function [ B ] = lila_ormqrf_v05( m, n, k, i, A, ia, ja, lda, B, ib, jb, ldb, T, mt )
%
      vb= mt - mod(i-1,mt) ;
      if ( vb > k ), vb = k; end
%
      ialo = ia+i-1;
      iahi = ia+m-1;
      jalo = ja+i-1;
      jahi = jalo+vb-1;
%
      iblo = ib+i-1;
      ibhi = ib+m-1;
      jblo = jb;
      jbhi = jb+n-1;
%
      itlo = mod(i-1,mt)+1;
      ithi = itlo+vb-1;
      jtlo = i;
      jthi = jtlo+vb-1;
%
      ml = m-i+1;
%
      not_done = 1;
%
      jj = 1;
%
      while ( not_done == 1 )

         V = tril(A(ialo:iahi,jalo:jahi), -1) + eye(ml,vb);

         H = (eye(ml,ml) - V * ( T(itlo:ithi,jtlo:jthi)' * V' ) );

         B(iblo:ibhi,jblo:jbhi) = H*B(iblo:ibhi,jblo:jbhi);
       
         if ( jj + vb - 1 == k ) 

            not_done = 0; 

         else if ( jj + vb - 1 > k ) 

             fprintf('defensive programming, we should never have been there, abort\n'); return;

         else

            ml = ml-vb;

            ialo = ialo+vb;
            iblo = iblo+vb;
            jalo = jalo+vb;
            itlo = 1;
            jtlo = jtlo+vb;

            jj = jj + vb;
 
            if ( ( jj+mt-1 ) <= k ), vb = mt; else vb = k-jj+1; end;

            jahi = jahi+vb;
            ithi = vb;
            jthi = jthi+vb;

         end
%
      end
%
   end
