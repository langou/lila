function [ A ] = lila_orgqr_v01( m, n, A, ia, ja, lda, T, it, jt, ldt )
%
      ialo = ia;
      iahi = ia+m-1;
      jalo = ja;
      jahi = ja+n-1;
%
      itlo = it;
      ithi = it+n-1;
      jtlo = jt;
      jthi = jt+n-1;
%
      QQ = zeros( m, n );
      QQ(1:n,1:n) = eye( n, n );
%
      V = tril( A( ialo:iahi, jalo:jahi ), -1 ) + eye( m, n );
      QQ = QQ - V * ( T(itlo:ithi,jtlo:jthi) * ( V' * QQ ) );
%
      A( ialo:iahi, jalo:jahi ) = QQ;
%
end
