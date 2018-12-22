%
   clear
%
   nb = [ 13, 10 ];
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
   ib = [ 4, 6, 3, 10 ];
   ib_block = ib_block + 1;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Construct Q
%
   QQ2 = Q(1:m,1:ihi(k-1));
   QQ2(1:m,ilo(k):ihi(k)) = zeros(m,nb(k));
   QQ2(ilo(k):ihi(k),ilo(k):ihi(k)) = eye(nb(k),nb(k));
%
   jlo = ihi(k)-ib(ib_block)+1;
   jhi = ihi(k);
%
%   This computes the full block of the "new" portion of Q
%   to allow us to focus on applying the little T's on the zeros of the new Q
%
%        The first time we apply the full new block QQ2 is the identity, so the calculation is not needed
%        Though if we don't do the full block in it's entirety this will need to be thought through
%
%         BBB(jlo:jhi,ilo(k):ihi(k)) = (tril(A(jlo:jhi,jlo:jhi),-1)+eye(ib(ib_block),ib(ib_block)))' * QQ2(jlo:jhi,ilo(k):ihi(k));
         BBB(jlo:jhi,ilo(k):ihi(k)) = (tril(A(jlo:jhi,jlo:jhi),-1)+eye(ib(ib_block),ib(ib_block)))';
%
         BBB(jlo:jhi,ilo(k):ihi(k)) = BBB(jlo:jhi,ilo(k):ihi(k)) + A(jhi+1:m,jlo:jhi)' * QQ2(jhi+1:m,ilo(k):ihi(k));
%
         BBB(jlo:jhi,ilo(k):ihi(k)) = T(jlo:jhi,jlo:jhi) * BBB(jlo:jhi,ilo(k):ihi(k));
%
         QQ2(jlo:jhi,ilo(k):ihi(k)) = QQ2(jlo:jhi,ilo(k):ihi(k)) - (tril(A(jlo:jhi,jlo:jhi),-1)+eye(ib(ib_block),ib(ib_block))) * BBB(jlo:jhi,ilo(k):ihi(k));
         QQ2(jhi+1:m,ilo(k):ihi(k)) = QQ2(jhi+1:m,ilo(k):ihi(k)) - A(jhi+1:m,jlo:jhi) * BBB(jlo:jhi,ilo(k):ihi(k));
%
         jlo = jlo - ib(ib_block-1);
         jhi = jhi - ib(ib_block);
%
%      for j = ib_block-1:-1:1,
%
%%%%%%%%%    Second step through ib to construct Q
%
%       This portion just continually put zeros in a part of BBB so I removed them
%
%         BBB(jlo:jhi,ilo(k):ihi(k)) = (tril(A(jlo:jhi,jlo:jhi),-1)+eye(ib(3),ib(3)))' * QQ2(jlo:jhi,ilo(k):ihi(k))
         BBB(jlo:jhi,ilo(k):ihi(k)) = zeros(ib(3),nb(k));

%
%       This doesn't need to addition of BBB anylonger since it's zero
%       A and Q are full and cannot make any saves
%
%         BBB(jlo:jhi,ilo(k):ihi(k)) = BBB(jlo:jhi,ilo(k):ihi(k)) + A(jhi+1:m,jlo:jhi)' * QQ2(jhi+1:m,ilo(k):ihi(k))
         BBB(jlo:jhi,ilo(k):ihi(k)) = A(jhi+1:m,jlo:jhi)' * QQ2(jhi+1:m,ilo(k):ihi(k));

%
%       Nothing to be saved flop-wise here 
%
         BBB(jlo:jhi,ilo(k):ihi(k)) = T(jlo:jhi,jlo:jhi) * BBB(jlo:jhi,ilo(k):ihi(k));
%
         QQ2(jlo:jhi,ilo(k):ihi(k)) = QQ2(jlo:jhi,ilo(k):ihi(k)) - (tril(A(jlo:jhi,jlo:jhi),-1)+eye(ib(3),ib(3))) * BBB(jlo:jhi,ilo(k):ihi(k));
         QQ2(jhi+1:m,ilo(k):ihi(k)) = QQ2(jhi+1:m,ilo(k):ihi(k)) - A(jhi+1:m,jlo:jhi) * BBB(jlo:jhi,ilo(k):ihi(k));
%
         jlo = jlo - ib(2); 
         jhi = jhi - ib(3);

%%%%%%%%%%     Third step through ib to construct Q

%         BBB(jlo:jhi,ilo(k):ihi(k)) = (tril(A(jlo:jhi,jlo:jhi),-1)+eye(ib(2),ib(2)))' * QQ2(jlo:jhi,ilo(k):ihi(k))
         BBB(jlo:jhi,ilo(k):ihi(k)) = zeros(ib(2),nb(k));

%         BBB(jlo:jhi,ilo(k):ihi(k)) = BBB(jlo:jhi,ilo(k):ihi(k)) + A(jhi+1:m,jlo:jhi)' * QQ2(jhi+1:m,ilo(k):ihi(k))
         BBB(jlo:jhi,ilo(k):ihi(k)) = A(jhi+1:m,jlo:jhi)' * QQ2(jhi+1:m,ilo(k):ihi(k));
%
%
%
         BBB(jlo:jhi,ilo(k):ihi(k)) = T(jlo:jhi,jlo:jhi) * BBB(jlo:jhi,ilo(k):ihi(k));
%
         QQ2(jlo:jhi,ilo(k):ihi(k)) = QQ2(jlo:jhi,ilo(k):ihi(k)) - (tril(A(jlo:jhi,jlo:jhi),-1)+eye(ib(2),ib(2))) * BBB(jlo:jhi,ilo(k):ihi(k));
         QQ2(jhi+1:m,ilo(k):ihi(k)) = QQ2(jhi+1:m,ilo(k):ihi(k)) - A(jhi+1:m,jlo:jhi) * BBB(jlo:jhi,ilo(k):ihi(k));
%

         jlo = jlo - ib(1); 
         jhi = jhi - ib(2); 
%
%%%%%%%%%

%         BBB(jlo:jhi,ilo(k):ihi(k)) = (tril(A(jlo:jhi,jlo:jhi),-1)+eye(ib(1),ib(1)))' * QQ2(jlo:jhi,ilo(k):ihi(k))
         BBB(jlo:jhi,ilo(k):ihi(k)) = zeros(ib(1),nb(k));

%         BBB(jlo:jhi,ilo(k):ihi(k)) = BBB(jlo:jhi,ilo(k):ihi(k)) + A(jhi+1:m,jlo:jhi)' * QQ2(jhi+1:m,ilo(k):ihi(k))
         BBB(jlo:jhi,ilo(k):ihi(k)) = A(jhi+1:m,jlo:jhi)' * QQ2(jhi+1:m,ilo(k):ihi(k));
%
%
%
         BBB(jlo:jhi,ilo(k):ihi(k)) = T(jlo:jhi,jlo:jhi) * BBB(jlo:jhi,ilo(k):ihi(k));
%
         QQ2(jlo:jhi,ilo(k):ihi(k)) = QQ2(jlo:jhi,ilo(k):ihi(k)) - (tril(A(jlo:jhi,jlo:jhi),-1)+eye(ib(1),ib(1))) * BBB(jlo:jhi,ilo(k):ihi(k));
         QQ2(jhi+1:m,ilo(k):ihi(k)) = QQ2(jhi+1:m,ilo(k):ihi(k)) - A(jhi+1:m,jlo:jhi) * BBB(jlo:jhi,ilo(k):ihi(k));
%
%%%%%%%%%

%         if (j > 1 ) jlo = jlo - ib(j-1); end;
%         if (j > 1 ) jhi = jhi - ib(j); end;
%

%       end

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
