%
   function [ A ] = lila_orgqr_v05( m, n, i, A, ia, ja, lda, T, it, jt, ldt )
%
     ml = m - i + 1 ;
%
     ialo = ia+i-1;
     iahi = m; % or we can write it as ( iahi = ia+i-1+ml-1; )
%
     jalo = ja+i-1;
     jahi = ja+i-1+n-1;
%
     itlo = it+i-1;
     ithi = it+i-1+n-1;
     jtlo = jt+i-1;
     jthi = jt+i-1+n-1;
%
     QQ = zeros( m, n );
     QQ(ialo:ialo+n-1,1:n) = eye( n,n );
%
     V = tril( A( ialo:iahi, jalo:jahi ), -1 ) + eye( ml, n );
T(itlo:ithi,jtlo:jthi) = larft( V );
     QQ(ialo:iahi,1:n) = QQ(ialo:iahi,1:n) - V * ( T(itlo:ithi,jtlo:jthi) * ( V' * QQ(ialo:iahi,1:n) ) );
%
     A( ialo:iahi, jalo:jahi ) = QQ(ialo:iahi,1:n);
%
     A( ialo:iahi, jalo:jahi ) = QQ( ialo:iahi, 1:n );
%
   end
