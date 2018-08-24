function [ A ] = lila_orgqr_v01( m, n, A, ia, ja, lda, T, AA )

      ilo = ia+n-1;
      ihi = ia+m-1;
      jlo = ja+n-1;
      jhi = ja+n-1;

      QQ = zeros( m, n );
      QQ(1:n,1:n) = eye( n, n );
%      for j = n:-1:1,
%         [ QQ(j:m,j:n) ] = larfL( A(ilo:ihi,jlo:jhi), QQ(j:m,j:n) );
%         ilo = ilo - 1;
%         jlo = jlo - 1;
%      end

%      A( ia:ihi, ja:jhi ) = QQ;


%
%  This is my comments for where I am stuck
%
%  So we didn't feed A to pull V - I've added that
%  Since we're working on the smaller Q the dimension
%  of H and QQ don't match. Do I need to only take a 
%  smaller H? Or should we go straight to projecting 
%  onto the larger space and construct the Q block
%  with [ zeros; eye; zeros ]?

      [mm,nn] = size(A);

      V = tril(AA(1:mm,1:ja), -1) + eye(mm,ja);

%  For T below - I know we need it to be square but the index I give
%  makes sense. How often would we run into a case where ia and ja are
%  not the same point? Since we're moving in blocks which end up being
%  square I don't think we'll need to worry about it here.

      H = (eye(mm,mm) - V * ( T(1:ia,1:ja) * V' ) );

size(H), size(QQ)

%      QQ(1:m,ja:jhi) = H*QQ(1:m,ja:jhi);
      QQ(1:m,ja:jhi) = H(1:m,1:m)*QQ(1:m,ja:jhi);



end
