%
   function [ Q ] = lila_ormqrbz_v05_w03( m, n, k, i, jb, mt, A, T, Q )
%
   iblo = i;
%
   jblo = jb;
   jbhi = jb+n-1;
%
   Q(iblo:i+k-1,jblo:jbhi) = zeros( size( Q(i:i+k-1,jblo:jbhi) ));
%
%
%
      vb= mt - mod(i-1,mt) ;
      if ( vb > k ), vb = k; end
%
      ialo = i+k-vb;
      iahi = m;
   
      jalo = i+k-vb;
      jahi = i+k-1;
  
      iblo = i+k-vb;
      ibhi = m;
   
      jblo = jb;
      jbhi = jb+n-1;
%
      itlo = mod(i-1,mt)+1;
      ithi = itlo-vb-1;
      jtlo = i;
      jthi = jtlo-vb-1;
%
      not_done = 1;
%
      jj = 1;
%
      while ( not_done == 1 )
%
         ml = iahi-ialo+1;
%
         V = tril( A( ialo:iahi,jalo:jahi ), -1 ) + eye( size( A( ialo:iahi,jalo:jahi ) ) );
% 
%         H = ( eye( ml,ml ) - V * ( T( itlo:ithi,jtlo:jthi ) * V' ) );
         H = ( eye( ml,ml ) - V * ( lapack_larft(V) * V' ) );
%
         Q( iblo:ibhi,jblo:jbhi ) = H * Q( iblo:ibhi,jblo:jbhi );
% 
         if ( jj + vb - 1 == k ) 

            not_done = 0; 

         else if ( jj + vb - 1 > k ) 

             fprintf('defensive programming, we should never have been there, abort\n'); return;

         else

            ml = ml+vb;

            ialo = ialo-vb;
            iblo = iblo-vb;
            jalo = jalo-vb;
            itlo = itlo-vb;
            jtlo = jtlo-vb;

            jj = jj + vb;
 
            if ( ( jj+mt-1 ) <= k ), vb = mt; else vb = k-jj+1; end;

            jahi = jahi-vb;
            ithi = vb;
            jthi = jthi-vb;

         end
%
      end
%
   end
