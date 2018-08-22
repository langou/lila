%
   clear
%
   nb = [ 19, 17 ];
   ib = [ 4, 9, 3, 3  ];
   m = 51;
%
   n = sum(nb);
   nb_block = size(nb,2);
   ib_block = size(ib,2);
%
   log10KA = 10;
%
   U = randn(m,n); [U,~]=qr(U,0);
   V = randn(n,n); [V,~]=qr(V,0);
   S = diag( 10.^( linspace( 0, -log10KA, n ) ) );
   A = U * S * V';
   clear U S V;
%
   ilo = zeros(nb_block,1);
   ilo(1) = 1;
   for i = 2:nb_block,
      ilo(i) = ilo(i-1) + nb(i-1);
   end
%
   ihi = zeros(nb_block,1);
   ihi = nb(1);
   for i = 2:nb_block,
      ihi(i) = ihi(i-1) + nb(i);
   end
%
   As = A;
   nrmA = norm(A,'fro');
   R = zeros(n);
   Q = zeros(m,n);
%
%%%%%%%%%%
%%%%%%%%%% Step 1 - can't be done in loop
%%%%%%%%%%
%
   [ A(ilo(1):m,ilo(1):ihi(1)) ] = geqr2( A(ilo(1):m,ilo(1):ihi(1)) );
%
   [ T(ilo(1):ihi(1),ilo(1):ihi(1)) ] = larft( A(ilo(1):m,ilo(1):ihi(1)) );
%
   Q(ilo(1):m,ilo(1):ihi(1)) = eye(m,nb(1));
   Q(ilo(1):m,ilo(1):ihi(1)) = Q(ilo(1):m,ilo(1):ihi(1)) - (tril(A(ilo(1):m,ilo(1):ihi(1)),-1)+eye(m,nb(1)))*(T(ilo(1):ihi(1),ilo(1):ihi(1))*((tril(A(ilo(1):m,ilo(1):ihi(1)),-1)+eye(m,nb(1)))'*Q(ilo(1):m,ilo(1):ihi(1))));
%
   ib = [ 4, 9, 3, 3  ];
   T( 1:ib(1), ib(1)+1:nb(1) ) = zeros( ib(1), nb(1)-ib(1) );
   T( ib(1)+1:ib(1)+ib(2), ib(1)+ib(2)+1:nb(1) ) = zeros( ib(2), nb(1)-ib(1)-ib(2) );
   T( ib(1)+ib(2)+1:ib(1)+ib(2)+ib(3), ib(1)+ib(2)+ib(3)+1:nb(1) ) = zeros( ib(3), nb(1)-ib(1)-ib(2)-ib(3) );


%%%%%%%%% Step 2:whatever 
%
   for k = 2;
%
   jlo = 1;
   jhi = ib(1);
   for j = 1:ib_block,
%
   BBB(jlo:jhi,ilo(k):ihi(k)) = A(jlo:jhi,ilo(k):ihi(k));
   BBB(jlo:jhi,ilo(k):ihi(k)) = (tril(A(jlo:jhi,jlo:jhi),-1)+eye(ib(j),ib(j)))'*BBB(jlo:jhi,ilo(k):ihi(k));
   BBB(jlo:jhi,ilo(k):ihi(k)) = BBB(jlo:jhi,ilo(k):ihi(k)) + A(jhi+1:m,jlo:jhi)'*A(jhi+1:m,ilo(k):ihi(k));
%
   BBB(jlo:jhi,ilo(k):ihi(k)) = triu( T(jlo:jhi,jlo:jhi) )'*BBB(jlo:jhi,ilo(k):ihi(k));
%
   A(jhi+1:m,ilo(k):ihi(k)) = A(jhi+1:m,ilo(k):ihi(k)) - A(jhi+1:m,jlo:jhi) * BBB(jlo:jhi,ilo(k):ihi(k));
   BBB(jlo:jhi,ilo(k):ihi(k)) = - (tril(A(jlo:jhi,jlo:jhi),-1)+eye(ib(j),ib(j))) * BBB(jlo:jhi,ilo(k):ihi(k));
%  [ BBB ] = blas_trmm ( 'L', 'L', 'N', 'U', ib(j), nb(k), -1.0e+00, A, jlo, jlo, -1, BBB, jlo, ilo(k), -1 );
   A(jlo:jhi,ilo(k):ihi(k)) = A(jlo:jhi,ilo(k):ihi(k)) + BBB(jlo:jhi,ilo(k):ihi(k));
%
   if (j < ib_block ) jlo = jlo + ib(j); end;
   if (j < ib_block ) jhi = jhi + ib(j+1); end;
   end
%
%
   [ A(ilo(k):m,ilo(k):ihi(k)) ] = geqr2( A(ilo(k):m,ilo(k):ihi(k)) );
%
   [ T(ilo(k):ihi(k),ilo(k):ihi(k)) ] = larft( A(ilo(k):m,ilo(k):ihi(k)) );
%
%  T(1:ihi(k-1),ilo(k):ihi(k)) = A(ilo(k):ihi(k),1:ihi(k-1))' * (tril(A(ilo(k):ihi(k),ilo(k):ihi(k)),-1)+eye(nb(k),nb(k))) ;
%  T(1:ihi(k-1),ilo(k):ihi(k)) = T(1:ihi(k-1),ilo(k):ihi(k)) + A(ihi(k)+1:m,1:ihi(k-1))' * A(ihi(k)+1:m,ilo(k):ihi(k))  ;
%  T(1:ihi(k-1),ilo(k):ihi(k)) = - ( T(1:ihi(k-1),1:ihi(k-1)) * T(1:ihi(k-1),ilo(k):ihi(k)) * T(ilo(k):ihi(k),ilo(k):ihi(k)) ) ;
%
%   ib = [ 4, 6, 3, 7, 7 ];
%   ib_block = ib_block + 2;
%  ib = [ 4, 6, 3, 4, 3, 7 ];
%  ib_block = ib_block + 3;
%

% The first 9 (from the right) is the kill point
% needs to be broken into [ 3, 6 ]
%
%   ib = [ 4, 9, 3, 3, 10, 5, 2 ];
   ib = [ 4, 9, 3, 9, 4, 5, 2 ];
   ib_block = ib_block + 3;


   T( 1:ib(1), ib(1)+1:nb(1) ) = zeros( ib(1), nb(1)-ib(1) );
   T( ib(1)+1:ib(1)+ib(2), ib(1)+ib(2)+1:nb(1) ) = zeros( ib(2), nb(1)-ib(1)-ib(2) );
   T( ib(1)+ib(2)+1:ib(1)+ib(2)+ib(3), ib(1)+ib(2)+ib(3)+1:nb(1) ) = zeros( ib(3), nb(1)-ib(1)-ib(2)-ib(3) );
   T( ib(1)+ib(2)+ib(3)+1:ib(1)+ib(2)+ib(3)+ib(4), ib(1)+ib(2)+ib(3)+ib(4)+1:nb(1)+nb(2) ) = zeros( ib(4), nb(1)+nb(2)-ib(1)-ib(2)-ib(3)-ib(4) );
   T( ib(1)+ib(2)+ib(3)+ib(4)+1:ib(1)+ib(2)+ib(3)+ib(4)+ib(5), ib(1)+ib(2)+ib(3)+ib(4)+ib(5)+1:nb(1)+nb(2) ) = zeros( ib(5), nb(1)+nb(2)-ib(1)-ib(2)-ib(3)-ib(4)-ib(5) );
   T( ib(1)+ib(2)+ib(3)+ib(4)+ib(5)+1:ib(1)+ib(2)+ib(3)+ib(4)+ib(5)+ib(6), ib(1)+ib(2)+ib(3)+ib(4)+ib(5)+ib(6)+1:nb(1)+nb(2) ) = zeros( ib(6), nb(1)+nb(2)-ib(1)-ib(2)-ib(3)-ib(4)-ib(5)-ib(6) );





%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Construct Q
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Constructing workspace
%
   QQ2 = Q(1:m,1:ihi(k-1));
%
   jlo = ihi(k)-ib(ib_block)+1;
   jhi = ihi(k);
%
%   (upper) T
%
      BBB(jlo:jhi,jlo:jhi) = (tril(A(jlo:jhi,jlo:jhi),-1)+eye(ib(ib_block),ib(ib_block)))';
%
%   (upper) T * (upper) T
%
      BBB(jlo:jhi,jlo:jhi) = triu( T(jlo:jhi,jlo:jhi) ) * triu( BBB(jlo:jhi,jlo:jhi) );
%
%    I - ( (lower) T * (upper) T )
%
      QQ2(jlo:jhi,jlo:jhi) = - (tril(A(jlo:jhi,jlo:jhi),-1)+eye(ib(ib_block),ib(ib_block))) * triu( BBB(jlo:jhi,jlo:jhi) );
      for i = jlo:jhi,
         QQ2(i,i) = 1.0e+00 + QQ2(i,i);
      end
%
%    Rect * (upper) T
%
      QQ2(jhi+1:m,jlo:jhi) = - A(jhi+1:m,jlo:jhi) * triu( BBB(jlo:jhi,jlo:jhi) );
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
      V=zeros(m,m);
      V(1:m,1:n) =  tril(A,-1)+eye(m,n);
%
      klo = ilo(k);
      khi = jhi - ib(ib_block);
%
      QQ2(klo:khi,klo:khi) = eye(size(QQ2(klo:khi,klo:khi)));
      QQ2(klo:khi,khi+1:ihi(k)) = zeros(size(QQ2(klo:khi,khi+1:ihi(k))));
      QQ2(khi+1:m,klo:khi) = zeros(size(QQ2(khi+1:m,klo:khi)));
      QQ2(khi+1:m,khi+1:ihi(k)) = QQ2(khi+1:m,khi+1:ihi(k));
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
      for j = ib_block-1:-1:ib_block-2,
%
%
      jlo = jlo - ib(j);
      jhi = jhi - ib(j+1);
%
      BBB(jlo:jhi,jlo:ihi(k)) = V(jlo:m,jlo:jhi)' * QQ2(jlo:m,jlo:ihi(k));
      BBB(jlo:jhi,ilo(k):ihi(k)) = T(jlo:jhi,jlo:jhi) * BBB(jlo:jhi,ilo(k):ihi(k));
      QQ2 ( jlo:m, ilo(k):ihi(k)) = QQ2 ( jlo:m, ilo(k):ihi(k)) - V(jlo:m,jlo:jhi) * BBB(jlo:jhi,ilo(k):ihi(k));
%
      end
%
%
      j = ib_block-3;
      ib1 = 6;
      ib2 = ib(j) - ib1;
      jlo = jlo - ib1;
      jhi = jhi - ib(j+1);
%
      BBB(jlo:jhi,jlo:ihi(k)) = V(jlo:m,jlo:jhi)' * QQ2(jlo:m,jlo:ihi(k));
      BBB(jlo:jhi,ilo(k):ihi(k)) = T(jlo:jhi,jlo:jhi) * BBB(jlo:jhi,ilo(k):ihi(k));
      QQ2 ( jlo:m, ilo(k):ihi(k)) = QQ2 ( jlo:m, ilo(k):ihi(k)) - V(jlo:m,jlo:jhi) * BBB(jlo:jhi,ilo(k):ihi(k));
%
      jlo = jlo - ib2;
      jhi = jhi - ib1;

      BBB(jlo:jhi,ilo(k):ihi(k)) = A(jhi+1:m,jlo:jhi)' * QQ2(jhi+1:m,ilo(k):ihi(k));
%
      BBB(jlo:jhi,ilo(k):ihi(k)) = T(jlo:jhi,jlo:jhi) * BBB(jlo:jhi,ilo(k):ihi(k));
%
      QQ2(jlo:jhi,ilo(k):ihi(k)) = QQ2(jlo:jhi,ilo(k):ihi(k)) - (tril(A(jlo:jhi,jlo:jhi),-1)+eye(ib2,ib2)) * BBB(jlo:jhi,ilo(k):ihi(k));
      QQ2(jhi+1:m,ilo(k):ihi(k)) = QQ2(jhi+1:m,ilo(k):ihi(k)) - A(jhi+1:m,jlo:jhi) * BBB(jlo:jhi,ilo(k):ihi(k));
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%    This portion constructs Q beyond the new block we've done QR on
%
      for j = ib_block-4:-1:1,
%
      jlo = jlo - ib(j);

      if ( j == ib_block-4 ) jhi = jhi - ib2; end;
      if ( j <  ib_block-4 ) jhi = jhi - ib(j+1); end;
%
         BBB(jlo:jhi,ilo(k):ihi(k)) = A(jhi+1:m,jlo:jhi)' * QQ2(jhi+1:m,ilo(k):ihi(k));
%
         BBB(jlo:jhi,ilo(k):ihi(k)) = T(jlo:jhi,jlo:jhi) * BBB(jlo:jhi,ilo(k):ihi(k));
%
         QQ2(jlo:jhi,ilo(k):ihi(k)) = QQ2(jlo:jhi,ilo(k):ihi(k)) - (tril(A(jlo:jhi,jlo:jhi),-1)+eye(ib(j),ib(j))) * BBB(jlo:jhi,ilo(k):ihi(k));
         QQ2(jhi+1:m,ilo(k):ihi(k)) = QQ2(jhi+1:m,ilo(k):ihi(k)) - A(jhi+1:m,jlo:jhi) * BBB(jlo:jhi,ilo(k):ihi(k));
%
%
       end
%
Q = QQ2;
%
   end
%
%
%
   R = triu( A(1:n,1:n) );
%
   fprintf('|| A - Q*R || / || A ||   = %6.1e\n',norm(As(1:m,1:n)-Q(1:m,1:n)*R(1:n,1:n),'fro')/nrmA);
   fprintf('|| I - Q''*Q ||            = %6.1e\n',norm(eye(n) - Q(1:m,1:n)'*Q(1:m,1:n),'fro'));
