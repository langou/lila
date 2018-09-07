%
   function [ B ] = lila_ormqrbz_v04( m, n, k, i, A, ia, ja, lda, B, ib, jb, ldb, T, it, jt, ldt, tb )
%
   Q = B;

   iblo = ib+i-1;

   jblo = jb;
   jbhi = jb+n-1;

   Q(iblo:i+k-1,jblo:jbhi) = zeros( size( Q(i:i+k-1,jblo:jbhi) ));

%%%%%%%%%%    This portion constructs Q column by column
%
%   ml = m - i - k + 2 ;
% 
%   ialo = ia+i-1+k-1;
%   iahi = m; % or we can write it as ( iahi = ia+i-1+ml-1; ) 
%  
%   jalo = ja+i-1+k-1;
%   jahi = ja+i-1+k-1;
% 
%   iblo = ib+i-1+k-1;
%   ibhi = m; % or we can write it as ( iahi = ia+i-1+ml-1; ) 
%  
%   jblo = jb;
%   jbhi = jb+n-1;
% 
%   itlo = it+i-1+k-1;
%   ithi = it+i-1+k-1;
%  
%   jtlo = jt+i-1+k-1;
%   jthi = jt+i-1+k-1;
%
%   for ii = i+k-1:-1:i,
%
%      V = tril(A(ialo:iahi,jalo:jahi), -1) + eye(size(A(ialo:iahi,jalo:jahi)));
%
%      H = (eye(ml,ml) - V * ( T(itlo:ithi,jtlo:jthi) * V' ) );
%
%      Q(iblo:ibhi,jblo:jbhi) = H * Q(iblo:ibhi,jblo:jbhi);
%
%      iblo = iblo-1;
%
%      ialo = ialo-1;
%      jalo = jalo-1;
%      jahi = jahi-1;
%
%      itlo = itlo-1;
%      ithi = ithi-1;
%
%      jtlo = jtlo-1;
%      jthi = jthi-1;
%
%      ml = ml+1;
%
%   end
%
%%%%%%%%%%    This portion constructs Q by a block of size - k   
%
%   ml = m - i + 1 ;
% 
%   ialo = ia+i-1;
%   iahi = m;  
%  
%   jalo = ja+i-1;
%   jahi = ja+i-1+k-1;
% 
%   iblo = ib+i-1;
%   ibhi = m;  
%  
%   jblo = jb;
%   jbhi = jb+n-1;
% 
%   itlo = it+i-1;
%   ithi = it+i-1+k-1;
%  
%   jtlo = jt+i-1;
%   jthi = jt+i-1+k-1;
%
%   V = tril( A( ialo:iahi,jalo:jahi ), -1 ) + eye( size( A( ialo:iahi,jalo:jahi ) ) );
%
%   H = ( eye( ml,ml ) - V * ( T( itlo:ithi,jtlo:jthi ) * V' ) );
%
%   Q( iblo:ibhi,jblo:jbhi ) = H * Q( iblo:ibhi,jblo:jbhi );
%
%%%%%%%%%%    This portion will construct Q based on tb    
%
   ll = 1;
   jj = 1;
   while( jj + tb( ll ) < i + n )
     jj = jj + tb( ll );
     ll = ll + 1;
   end
   vb = i + n - jj - 2;  %Why does minus two work?
%
%
   ialo = ia+i-1+k-1-vb+1;
   iahi = m; 
% 
   jalo = ja+i-1+k-1-vb+1;
   jahi = ja+i-1+k-1;
%
   itlo = it+i-1+k-1-vb+1;
   ithi = it+i-1+k-1;
% 
   jtlo = jt+i-1+k-1-vb+1;
   jthi = jt+i-1+k-1;
%
   iblo = ib+i-1+k-1-vb+1;
   ibhi = m;
%
   jblo = jb;
   jbhi = jb+n-1;
%
   not_done = 1;
%
   while ( not_done == 1 )
%
      ml = iahi-ialo+1;
%
      V = tril( A( ialo:iahi,jalo:jahi ), -1 ) + eye( size( A( ialo:iahi,jalo:jahi ) ) );
%
      H = ( eye( ml,ml ) - V * ( T( itlo:ithi,jtlo:jthi ) * V' ) );
%
      Q( iblo:ibhi,jblo:jbhi ) = H * Q( iblo:ibhi,jblo:jbhi );
% 
      jahi = jahi-vb;     
      ithi = ithi-vb; 
      jthi = jthi-vb;    
%
      ll = ll - 1;
      if( ll > 0 )
      if ( i <= ( ialo-tb(ll) ) ), vb = tb(ll); else vb = ialo - i; not_done = 0; end;
      else 
         not_done = 0;
      end
%
      ialo = ialo-vb;
      jalo = jalo-vb;         
      itlo = itlo-vb;          
      jtlo = jtlo-vb;     
      iblo = iblo-vb;    
%
   end
%
   B = Q;
%
   end
