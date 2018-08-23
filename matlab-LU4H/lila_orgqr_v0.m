function [ A ] = lila_orgqr_v0( m, n, A, ia, ja, lda )


      QQ = zeros( ml, nn );
      QQ(1:nn,1:nn) = eye( nn, nn );
      for j = nn:-1:1,
         [ QQ(j:ml,j:nn) ] = larfL( AA(j:ml,j), QQ(j:ml,j:nn) );
      end

% copy QQ back in A AT THE RIGHT SPOT
% A( ilo:ihi, jlo:jhi ) = QQ:


