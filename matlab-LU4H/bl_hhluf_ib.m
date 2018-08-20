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
% 1
   [ A(ilo(1):m,ilo(1):ihi(1)) ] = geqr2( A(ilo(1):m,ilo(1):ihi(1)) );
%
   [ T(ilo(1):ihi(1),ilo(1):ihi(1)) ] = larft( A(ilo(1):m,ilo(1):ihi(1)) );
%
   Q(ilo(1):m,ilo(1):ihi(1)) = eye(m,nb(1));
   Q(ilo(1):m,ilo(1):ihi(1)) = Q(ilo(1):m,ilo(1):ihi(1)) - (tril(A(ilo(1):m,ilo(1):ihi(1)),-1)+eye(m,nb(1)))*(T(ilo(1):ihi(1),ilo(1):ihi(1))*((tril(A(ilo(1):m,ilo(1):ihi(1)),-1)+eye(m,nb(1)))'*Q(ilo(1):m,ilo(1):ihi(1))));
%
%
   for k = 2;
%
%  AA1 = A;
%  AA2 = A;
%  BBB = ((tril(A(1:ihi(k-1),1:ihi(k-1)),-1)+eye(ihi(k-1),ihi(k-1)))'*A(1:ihi(k-1),ilo(k):ihi(k)));
%  BBB = BBB + A(ilo(k):m,1:ihi(k-1))'*A(ilo(k):m,ilo(k):ihi(k));
%  BBB = triu( T(1:ihi(k-1),1:ihi(k-1)) )'*BBB;
%  A(1:ihi(k-1),ilo(k):ihi(k)) = A(1:ihi(k-1),ilo(k):ihi(k)) - (tril(A(1:ihi(k-1),1:ihi(k-1)),-1)+eye(ihi(k-1),ihi(k-1)))*BBB;
%  A(ilo(k):m,ilo(k):ihi(k)) = A(ilo(k):m,ilo(k):ihi(k)) - A(ilo(k):m,1:ihi(k-1))*BBB;
%  A1 = A;
%
%  V = (tril(AA1(1:m,1:ihi(k-1)),-1)+eye(m,ihi(k-1)));
%  BB1 = V'*AA1(1:m,ilo(k):ihi(k));
%  BB1 = triu( T(1:ihi(k-1),1:ihi(k-1)) )'*BB1;
%  AA1(1:m,ilo(k):ihi(k)) = AA1(1:m,ilo(k):ihi(k)) - V * BB1;
%
   jlo = 1;
   jhi = ib(1);
   for j = 1:ib_block,

   BBB(jlo:jhi,ilo(k):ihi(k)) = A(jlo:jhi,ilo(k):ihi(k));
   BBB(jlo:jhi,ilo(k):ihi(k)) = (tril(A(jlo:jhi,jlo:jhi),-1)+eye(ib(j),ib(j)))'*BBB(jlo:jhi,ilo(k):ihi(k));
   BBB(jlo:jhi,ilo(k):ihi(k)) = BBB(jlo:jhi,ilo(k):ihi(k)) + A(jhi+1:m,jlo:jhi)'*A(jhi+1:m,ilo(k):ihi(k));

   BBB(jlo:jhi,ilo(k):ihi(k)) = triu( T(jlo:jhi,jlo:jhi) )'*BBB(jlo:jhi,ilo(k):ihi(k));

   A(jhi+1:m,ilo(k):ihi(k)) = A(jhi+1:m,ilo(k):ihi(k)) - A(jhi+1:m,jlo:jhi) * BBB(jlo:jhi,ilo(k):ihi(k));
   BBB(jlo:jhi,ilo(k):ihi(k)) = - (tril(A(jlo:jhi,jlo:jhi),-1)+eye(ib(j),ib(j))) * BBB(jlo:jhi,ilo(k):ihi(k));
%  [ BBB ] = blas_trmm ( 'L', 'L', 'N', 'U', ib(j), nb(k), -1.0e+00, A, jlo, jlo, -1, BBB, jlo, ilo(k), -1 );
   A(jlo:jhi,ilo(k):ihi(k)) = A(jlo:jhi,ilo(k):ihi(k)) + BBB(jlo:jhi,ilo(k):ihi(k));

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
%  ib = [ 4, 6, 3, 3, 2, 5  ];
%  ib = [ 4, 6, 5, 2, 2, 4  ];
   ib = [ 4, 6, 3, 2, 2, 2, 4  ];
   ib_block = ib_block + 4;
%
   BBB(1:ihi(k-1),1:nb(k)) = triu( T(1:ihi(k-1),1:ihi(k-1)) )*[A(ilo(k):ihi(k),1:ihi(k-1))'];
   BBB(1:ihi(k-1),1:nb(k)) = BBB(1:ihi(k-1),1:nb(k)) + T(1:ihi(k-1),ilo(k):ihi(k))*( tril(A(ilo(k):ihi(k),ilo(k):ihi(k)),-1)+eye(nb(k)) )';
   BBB(ilo(k):ihi(k),1:nb(k)) = triu( T(ilo(k):ihi(k),ilo(k):ihi(k)) )*[tril(A(ilo(k):ihi(k),ilo(k):ihi(k)),-1)+eye(nb(k))]';
%
   Q(1:ihi(k-1),ilo(k):ihi(k)) = - (tril(A(1:ihi(k-1),1:ihi(k-1)),-1)+eye(ihi(k-1),ihi(k-1))) * BBB(1:ihi(k-1),1:nb(k));
%
   Q(ilo(k):ihi(k),ilo(k):ihi(k)) = ( tril( A(ilo(k):ihi(k),ilo(k):ihi(k)), -1 ) + eye(nb(k),nb(k) ) ) * triu( BBB(ilo(k):ihi(k),1:nb(k)) ) ;
   Q(ilo(k):ihi(k),ilo(k):ihi(k)) = eye(nb(k),nb(k)) - Q(ilo(k):ihi(k),ilo(k):ihi(k));
   Q(ilo(k):ihi(k),ilo(k):ihi(k)) = Q(ilo(k):ihi(k),ilo(k):ihi(k)) - A(ilo(k):ihi(k),1:ihi(k-1)) * BBB(1:ihi(k-1),1:nb(k));
%
   Q(ihi(k)+1:m,ilo(k):ihi(k)) = - A(ihi(k)+1:m,ilo(k):ihi(k)) * triu( BBB(ilo(k):ihi(k),1:nb(k)) );
   Q(ihi(k)+1:m,ilo(k):ihi(k)) = Q(ihi(k)+1:m,ilo(k):ihi(k)) - A(ihi(k)+1:m,1:ihi(k-1)) * BBB(1:ihi(k-1),1:nb(k)) ;
%
   QQ2 = Q(1:m,1:ihi(k-1));
   QQ2(1:m,ilo(k):ihi(k)) = zeros(m,nb(k));
   QQ2(ilo(k):ihi(k),ilo(k):ihi(k)) = eye(nb(k),nb(k));

   V = (tril(A(1:m,1:ihi(k)),-1)+eye(m,ihi(k)));

%  QQ2(1:m,ilo(k):ihi(k)) = ( eye(m) - V(1:m,1:ihi(k)) * T(1:ihi(k),1:ihi(k)) * V(1:m,1:ihi(k))' ) * QQ2(1:m,ilo(k):ihi(k));

%  for j = ihi(k):-1:1,
%  QQ2(1:m,ilo(k):ihi(k)) = ( eye(m) - V(1:m,j) * T(j,j) * V(1:m,j)' ) * QQ2(1:m,ilo(k):ihi(k));
%  end

%  jlo = ihi(k)-ib(ib_block)+1;
%  jhi = ihi(k);
%  for j = ib_block:-1:1,
%
%  QQ2(1:m,ilo(k):ihi(k)) = ( eye(m) - V(1:m,jlo:jhi) * T(jlo:jhi,jlo:jhi) * V(1:m,jlo:jhi)' ) * QQ2(1:m,ilo(k):ihi(k));
%
%  if (j > 1 ) jlo = jlo - ib(j-1); end;
%  if (j > 1 ) jhi = jhi - ib(j); end;
%
%  end




%
   QQ2 = Q(1:m,1:ihi(k-1));
   QQ2(1:m,ilo(k):ihi(k)) = zeros(m,nb(k));
   QQ2(ilo(k):ihi(k),ilo(k):ihi(k)) = eye(nb(k),nb(k));




QQ2(1:m,ilo(k):ihi(k))


   jlo = ihi(k)-ib(ib_block)+1;
   jhi = ihi(k);
   for j = ib_block:-1:1,

   V = (tril(A(1:m,1:ihi(k)),-1)+eye(m,ihi(k)));

   BBB(jlo:jhi,ilo(k):ihi(k)) = (tril(A(jlo:jhi,jlo:jhi),-1)+eye(ib(j),ib(j)))' * QQ2(jlo:jhi,ilo(k):ihi(k));
   BBB(jlo:jhi,ilo(k):ihi(k)) = BBB(jlo:jhi,ilo(k):ihi(k)) + A(jhi+1:m,jlo:jhi)' * QQ2(jhi+1:m,ilo(k):ihi(k));

   BBB(jlo:jhi,ilo(k):ihi(k)) = T(jlo:jhi,jlo:jhi) * BBB(jlo:jhi,ilo(k):ihi(k));

   QQ2(jlo:jhi,ilo(k):ihi(k)) = QQ2(jlo:jhi,ilo(k):ihi(k)) - (tril(A(jlo:jhi,jlo:jhi),-1)+eye(ib(j),ib(j))) * BBB(jlo:jhi,ilo(k):ihi(k));
   QQ2(jhi+1:m,ilo(k):ihi(k)) = QQ2(jhi+1:m,ilo(k):ihi(k)) - A(jhi+1:m,jlo:jhi) * BBB(jlo:jhi,ilo(k):ihi(k));

QQ2(1:m,ilo(k):ihi(k))


   if (j > 1 ) jlo = jlo - ib(j-1); end;
   if (j > 1 ) jhi = jhi - ib(j); end;
 
   end





   end
%
%
%
   R = triu( A(1:n,1:n) );
%
   fprintf('|| A - Q*R || / || A ||   = %6.1e\n',norm(As(1:m,1:n)-Q(1:m,1:n)*R(1:n,1:n),'fro')/nrmA);
   fprintf('|| I - Q''*Q ||            = %6.1e\n',norm(eye(n) - Q(1:m,1:n)'*Q(1:m,1:n),'fro'));
