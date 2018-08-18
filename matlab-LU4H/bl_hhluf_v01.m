%
   clear
%
   nb = [ 2, 3, 4, 7, 5 ];
   m = 40;
%
   n = sum(nb);
   nb_block = size(nb,2);
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
%  3
   for k = 2:5;
   BBB = ((tril(A(1:ihi(k-1),1:ihi(k-1)),-1)+eye(ihi(k-1),ihi(k-1)))'*A(1:ihi(k-1),ilo(k):ihi(k)));
   BBB = BBB + A(ilo(k):m,1:ihi(k-1))'*A(ilo(k):m,ilo(k):ihi(k));
   BBB = triu( T(1:ihi(k-1),1:ihi(k-1)) )'*BBB;
   A(1:ihi(k-1),ilo(k):ihi(k)) = A(1:ihi(k-1),ilo(k):ihi(k)) - (tril(A(1:ihi(k-1),1:ihi(k-1)),-1)+eye(ihi(k-1),ihi(k-1)))*BBB;
   A(ilo(k):m,ilo(k):ihi(k)) = A(ilo(k):m,ilo(k):ihi(k)) - A(ilo(k):m,1:ihi(k-1))*BBB;
%
   [ A(ilo(k):m,ilo(k):ihi(k)) ] = geqr2( A(ilo(k):m,ilo(k):ihi(k)) );
%
   [ T(ilo(k):ihi(k),ilo(k):ihi(k)) ] = larft( A(ilo(k):m,ilo(k):ihi(k)) );
%
   T(1:ihi(k-1),ilo(k):ihi(k)) = A(ilo(k):ihi(k),1:ihi(k-1))' * (tril(A(ilo(k):ihi(k),ilo(k):ihi(k)),-1)+eye(nb(k),nb(k))) ;
   T(1:ihi(k-1),ilo(k):ihi(k)) = T(1:ihi(k-1),ilo(k):ihi(k)) + A(ihi(k)+1:m,1:ihi(k-1))' * A(ihi(k)+1:m,ilo(k):ihi(k))  ;
   T(1:ihi(k-1),ilo(k):ihi(k)) = - ( T(1:ihi(k-1),1:ihi(k-1)) * T(1:ihi(k-1),ilo(k):ihi(k)) * T(ilo(k):ihi(k),ilo(k):ihi(k)) ) ;
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
   end
%
%
%
   R = triu( A(1:n,1:n) );
%
   fprintf('|| A - Q*R || / || A ||   = %6.1e\n',norm(As(1:m,1:n)-Q(1:m,1:n)*R(1:n,1:n),'fro')/nrmA);
   fprintf('|| I - Q''*Q ||            = %6.1e\n',norm(eye(n) - Q(1:m,1:n)'*Q(1:m,1:n),'fro'));
