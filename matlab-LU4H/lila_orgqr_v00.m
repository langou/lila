function [ A ] = lila_orgqr_v0( m, n, A, ia, ja, lda )

      ilo = ia+n-1;
      ihi = ia+m-1;
      jlo = ja+n-1;
      jhi = ja+n-1;

      QQ = zeros( m, n );
      QQ(1:n,1:n) = eye( n, n );
      for j = n:-1:1,
         [ QQ(j:m,j:n) ] = larfL( A(ilo:ihi,jlo:jhi), QQ(j:m,j:n) );
         ilo = ilo - 1;
         jlo = jlo - 1;
      end

      A( ia:ihi, ja:jhi ) = QQ;

end
