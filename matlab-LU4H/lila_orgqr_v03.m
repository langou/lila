%
 function [ A ] = lila_orgqr_v03( m, n, i, A, ia, ja, lda, T, it, jt, ldt, tb )
%
      ml = m - i + 1 ;
%
%%%% this is the whole block at once (no transpose T)
%
%     ialo = ia+i-1;
%     iahi = m; % or we can write it as ( iahi = ia+i-1+ml-1; )
%
%     jalo = ja+i-1;
%     jahi = ja+i-1+n-1;
%
%     itlo = it+i-1;
%     ithi = it+i-1+n-1;
%     jtlo = jt+i-1;
%     jthi = jt+i-1+n-1;
%
%     QQ = zeros( m, n );
%     QQ(ialo:ialo+n-1,1:n) = eye( n,n );
%
%     V = tril( A( ialo:iahi, jalo:jahi ), -1 ) + eye( ml, n );
%     QQ(ialo:iahi,1:n) = QQ(ialo:iahi,1:n) - V * ( T(itlo:ithi,jtlo:jthi) * ( V' * QQ(ialo:iahi,1:n) ) );
%
%     A( ialo:iahi, jalo:jahi ) = QQ(ialo:iahi,1:n);
%
%%%% this is one-by-one from the bottom
%
%   mll = ml;
%   nll = n;               
% 
%   ialo = ia+i-1+n-1;
%   iahi = m; % or we can write it as ( iahi = ia+i-1+ml-1; ) 
%  
%   jalo = ja+i-1+n-1;
%   jahi = ja+i-1+n-1;
% 
%   itlo = it+i-1+n-1;
%   ithi = it+i-1+n-1;
%  
%   jtlo = jt+i-1+n-1;
%   jthi = jt+i-1+n-1;
% 
%   QQQQ = zeros( m, n );
%   QQ = zeros( ml, n );
%
%   QQQQ(i:i+n-1,1:n) = eye( n,n );
%   QQ(1:n,1:n) = eye( n,n );
%
%   for ii = n:-1:1,
%
%      V = tril( A( ialo:iahi, jalo:jahi ), -1 ) + eye( size(ialo:iahi,2),1 );
%
%      QQ( ii:mll,ii:nll ) = QQ( ii:mll,ii:nll ) - V * ( T(itlo:ithi,jtlo:jthi ) * ( V' * QQ( ii:mll,ii:nll ) ) );
%      QQQQ( ii+i-1:m,ii:nll ) = QQQQ( ii+i-1:m,ii:nll ) - V * ( T(itlo:ithi,jtlo:jthi ) * ( V' * QQQQ( ii+i-1:m,ii:nll ) ) );   
% 
%      ialo = ialo-1;
% 
%      jalo = jalo-1;         
%      jahi = jahi-1;     
% 
%      itlo = itlo-1;          
%      ithi = ithi-1; 
% 
%      jtlo = jtlo-1;         
%      jthi = jthi-1;    
% 
%   end
% 
%   ialo = ia+i-1;
%   iahi = ia+i-1+ml-1; 
%   jalo = ja+i-1;
%   jahi = ja+i-1+n-1;
%
%   A( ialo:iahi, jalo:jahi ) = QQ;
%   A( ialo:iahi, jalo:jahi ) = QQQQ( ialo:iahi, 1:n);
%
%%%% this is by blocks respecting the tb structure at starting from the bottom
%
%  tb = [ 4, 3, 9, 5, 4, 5, 7 ];
%  nb = [ 10, 13 ];
%
%  ub = [ 4, 3, 3, 6, 5, 2, 2, 5, 7 ];
%                  ^  ^  ^
%  ll =            4  5  6

   ll = 1;
   jj = 1;
   while( jj + tb( ll ) < i + n )
     jj = jj + tb( ll );
     ll = ll + 1;
   end
   vb = i + n - jj;

%
   QQ = zeros( m,n );
   QQ(i:i+n-1,1:n) = eye( n,n );

   mll = ml;               
   nll = n;               

   ialo = ia+i-1+n-1-vb+1;
   iahi = ia+i-1+ml-1; % this is m, isn't it?
 
   jalo = ja+i-1+n-1-vb+1;
   jahi = ja+i-1+n-1;

   itlo = it+i-1+n-1-vb+1;
   ithi = it+i-1+n-1;
 
   jtlo = jt+i-1+n-1-vb+1;
   jthi = jt+i-1+n-1;

   not_done = 1;

   while ( not_done == 1 )


      V = tril( A( ialo:iahi, jalo:jahi ), -1 ) + eye( size( A( ialo:iahi, jalo:jahi ) ));
  
      QQ(ialo:iahi,1:n) = QQ(ialo:iahi,1:n) - V * ( T(itlo:ithi,jtlo:jthi) * ( V' * QQ(ialo:iahi,1:n) ) );
  
      jahi = jahi-vb;     
      ithi = ithi-vb; 
      jthi = jthi-vb;    

      ll = ll - 1;
      if( ll > 0 )
      if ( i <= ( ialo-tb(ll) ) ), vb = tb(ll); else vb = ialo - i; not_done = 0; end;
      else 
         not_done = 0;
      end


      ialo = ialo-vb;
      jalo = jalo-vb;         
      itlo = itlo-vb;          
      jtlo = jtlo-vb;         

   end

%%%%%

  
       V = tril( A( ialo:iahi, jalo:jahi ), -1 ) + eye( size( A( ialo:iahi, jalo:jahi ) ));
   
       QQ(ialo:iahi,1:n) = QQ(ialo:iahi,1:n) - V * ( T(itlo:ithi,jtlo:jthi) * ( V' * QQ(ialo:iahi,1:n) ) );
   





  ialo = ia+i-1;
  iahi = ia+i-1+ml-1; % this is m, isn't it?
  jalo = ja+i-1;
  jahi = ja+i-1+n-1;

  A( ialo:iahi, jalo:jahi ) = QQ( ialo:iahi, 1:n );







end
