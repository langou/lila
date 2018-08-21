%
   clear
%
   nb = [ 13, 14 ];
   ib = [ 4, 6, 3  ];
   m = 40;
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
   T(1:ihi(k-1),ilo(k):ihi(k)) = A(ilo(k):ihi(k),1:ihi(k-1))' * (tril(A(ilo(k):ihi(k),ilo(k):ihi(k)),-1)+eye(nb(k),nb(k))) ;
   T(1:ihi(k-1),ilo(k):ihi(k)) = T(1:ihi(k-1),ilo(k):ihi(k)) + A(ihi(k)+1:m,1:ihi(k-1))' * A(ihi(k)+1:m,ilo(k):ihi(k))  ;
   T(1:ihi(k-1),ilo(k):ihi(k)) = - ( T(1:ihi(k-1),1:ihi(k-1)) * T(1:ihi(k-1),ilo(k):ihi(k)) * T(ilo(k):ihi(k),ilo(k):ihi(k)) ) ;
%
%   ib = [ 4, 6, 3, 7, 7 ];
%   ib_block = ib_block + 2;
%  ib = [ 4, 6, 3, 4, 3, 7 ];
%  ib_block = ib_block + 3;
%
   ib = [ 3, 5, 2, 9, 3, 5 ];
   ib_block = ib_block + 3;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Construct Q
%
%      Breaking 10 into two
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
%
      V=zeros(m,m);
      V(1:m,1:n) =  tril(A,-1)+eye(m,n);
%


      klo = jlo - ib(ib_block-1)- ib(ib_block-2);
      khi = jhi - ib(ib_block);

      QQ2(klo:khi,klo:khi) = eye(size(QQ2(klo:khi,klo:khi)));
      QQ2(klo:khi,khi+1:ihi(k)) = zeros(size(QQ2(klo:khi,khi+1:ihi(k))));
      QQ2(khi+1:m,klo:khi) = zeros(size(QQ2(khi+1:m,klo:khi)));
      QQ2(khi+1:m,khi+1:ihi(k)) = QQ2(khi+1:m,khi+1:ihi(k));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      k_break = 4; % 4 is the block we want to break

      for j = ib_block-1:-1:k_break+1,

      jlo = jlo - ib(j);
      jhi = jhi - ib(j+1);

      BBB(jlo:jhi,jlo:ihi(k)) = V(jlo:m,jlo:jhi)' * QQ2(jlo:m,jlo:ihi(k));
      BBB(jlo:jhi,ilo(k):ihi(k)) = T(jlo:jhi,jlo:jhi) * BBB(jlo:jhi,ilo(k):ihi(k));
      QQ2 ( jlo:m, ilo(k):ihi(k)) = QQ2 ( jlo:m, ilo(k):ihi(k)) - V(jlo:m,jlo:jhi) * BBB(jlo:jhi,ilo(k):ihi(k));

      end

      i1 = 3; % 3 is how we break block k_break so we break 9=3+6
      i2 = ib(k_break) - i1;
 
      jlo = jlo - i2;
      jhi = jhi - ib(k_break+1);

      BBB(jlo:jhi,jlo:ihi(k)) = V(jlo:m,jlo:jhi)' * QQ2(jlo:m,jlo:ihi(k));
      BBB(jlo:jhi,ilo(k):ihi(k)) = T(jlo:jhi,jlo:jhi) * BBB(jlo:jhi,ilo(k):ihi(k));
      QQ2 ( jlo:m, ilo(k):ihi(k)) = QQ2 ( jlo:m, ilo(k):ihi(k)) - V(jlo:m,jlo:jhi) * BBB(jlo:jhi,ilo(k):ihi(k));

      jlo = jlo - i1;
      jhi = jhi - i2;

         BBB(jlo:jhi,ilo(k):ihi(k)) = A(jhi+1:m,jlo:jhi)' * QQ2(jhi+1:m,ilo(k):ihi(k));
%
         BBB(jlo:jhi,ilo(k):ihi(k)) = T(jlo:jhi,jlo:jhi) * BBB(jlo:jhi,ilo(k):ihi(k));
%
         QQ2(jlo:jhi,ilo(k):ihi(k)) = QQ2(jlo:jhi,ilo(k):ihi(k)) - (tril(A(jlo:jhi,jlo:jhi),-1)+eye(i1,i1)) * BBB(jlo:jhi,ilo(k):ihi(k));
         QQ2(jhi+1:m,ilo(k):ihi(k)) = QQ2(jhi+1:m,ilo(k):ihi(k)) - A(jhi+1:m,jlo:jhi) * BBB(jlo:jhi,ilo(k):ihi(k));
%


 


 


%




%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% NOW
%
%    (upper) T
%
%     BBB(jlo:jhi,jlo:jhi) = (tril(A(jlo:jhi,jlo:jhi),-1)+eye(ib(ib_block-2),ib(ib_block-2)))';
%
%%%%%%   Column had zeros based off the previous ib
%
%   Rect * Rect
%
%     BBB(jlo:jhi,ilo(k)+ib(ib_block-1):ihi(k)) = A(jhi+1:m,jlo:jhi)' * QQ2(jhi+1:m,ilo(k)+ib(ib_block-1):ihi(k));
%%%%%%
%
%   (upper) T * Rect
%
%     BBB(jlo:jhi,ilo(k):ihi(k)) = T(jlo:jhi,jlo:jhi) * BBB(jlo:jhi,ilo(k):ihi(k));
%
%   [ I, zeros ]  - (lower) T * [ (upper) T, Rect ]
%
%     QQ2(jlo:jhi,ilo(k):ihi(k)) = QQ2(jlo:jhi,ilo(k):ihi(k)) - (tril(A(jlo:jhi,jlo:jhi),-1)+eye(ib(ib_block-2),ib(ib_block-2))) * BBB(jlo:jhi,ilo(k):ihi(k));
%
%   [ zeros, Rect ]  - Rect * [ (upper) T, Rect ]
%
%     QQ2(jhi+1:m,ilo(k):ihi(k)) = QQ2(jhi+1:m,ilo(k):ihi(k)) - A(jhi+1:m,jlo:jhi) * BBB(jlo:jhi,ilo(k):ihi(k));
%
%     jlo = jlo - ib(ib_block-3);
%     jhi = jhi - ib(ib_block-2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% NOW
%
%    This portion constructs Q beyond the new block we've done QR on
%
      for j = k_break-1:-1:1,

      jlo = jlo - ib(j);

      if (j==k_break-1) jhi = jhi - i1; else jhi = jhi - ib(j+1); end


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
