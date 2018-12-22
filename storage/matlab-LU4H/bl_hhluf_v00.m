%
   clear
%
   nb = [ 2, 3, 4 ];
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
   ilo(1) = [1];
   for i = 2:nb_block,
      ilo(i) = ilo(i-1) + nb(i-1);
   end
%
   ihi = [nb(1)];
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
%  2
   A(ilo(1):m,ilo(2):ihi(2)) = A(ilo(1):m,ilo(2):ihi(2)) - (tril(A(ilo(1):m,ilo(1):ihi(1)),-1)+eye(m,nb(1)))*(T(ilo(1):ihi(1),ilo(1):ihi(1))'*((tril(A(ilo(1):m,ilo(1):ihi(1)),-1)+eye(m,nb(1)))'*A(ilo(1):m,ilo(2):ihi(2))));
%
   [ A(ilo(2):m,ilo(2):ihi(2)) ] = geqr2( A(ilo(2):m,ilo(2):ihi(2)) );
%
   [ T(ilo(2):ihi(2),ilo(2):ihi(2)) ] = larft( A(ilo(2):m,ilo(2):ihi(2)) );
%
   Q(ilo(2):m,ilo(2):ihi(2)) = eye(m-nb(1),nb(2));
%
   Q(ilo(2):m,ilo(2):ihi(2)) = Q(ilo(2):m,ilo(2):ihi(2)) - (tril(A(ilo(2):m,ilo(2):ihi(2)),-1)+eye(m-nb(1),nb(2)))*(T(ilo(2):ihi(2),ilo(2):ihi(2))*((tril(A(ilo(2):m,ilo(2):ihi(2)),-1)+eye(m-nb(1),nb(2)))'*Q(ilo(2):m,ilo(2):ihi(2))));
%
   Q(ilo(1):m,ilo(2):ihi(2)) = Q(ilo(1):m,ilo(2):ihi(2)) - (tril(A(ilo(1):m,ilo(1):ihi(1)),-1)+eye(m,nb(1)))*(T(ilo(1):ihi(1),ilo(1):ihi(1))*((tril(A(ilo(1):m,ilo(1):ihi(1)),-1)+eye(m,nb(1)))'* ...
           [ zeros(nb(1),nb(2)); Q(ilo(2):m,ilo(2):ihi(2)) ] ));
%
   T(ilo(1):ihi(1),ilo(2):ihi(2)) = - ( T(ilo(1):ihi(1),ilo(1):ihi(1)) * ( A(ilo(2):m,ilo(1):ihi(1))' * (tril(A(ilo(2):m,ilo(2):ihi(2)),-1)+eye(m-nb(1),nb(2))) ) * T(ilo(2):ihi(2),ilo(2):ihi(2)) ) ;
%
   Q(1:m,ilo(2):ihi(2)) = zeros(m,nb(2));
   Q(ilo(2):m,ilo(2):ihi(2)) = eye(m-nb(1),nb(2));
%
   Q(ilo(2):m,ilo(2):ihi(2)) = Q(ilo(2):m,ilo(2):ihi(2)) - (tril(A(ilo(2):m,ilo(2):ihi(2)),-1)+eye(m-nb(1),nb(2)))*(T(ilo(2):ihi(2),ilo(2):ihi(2))*((tril(A(ilo(2):m,ilo(2):ihi(2)),-1)+eye(m-nb(1),nb(2)))'*Q(ilo(2):m,ilo(2):ihi(2))));
%
   Q(ilo(1):m,ilo(2):ihi(2)) = Q(ilo(1):m,ilo(2):ihi(2)) - (tril(A(ilo(1):m,ilo(1):ihi(1)),-1)+eye(m,nb(1)))*(T(ilo(1):ihi(1),ilo(1):ihi(1))*((tril(A(ilo(1):m,ilo(1):ihi(1)),-1)+eye(m,nb(1)))'* ...
           [ zeros(nb(1),nb(2)); Q(ilo(2):m,ilo(2):ihi(2)) ] ));
%
%  3
   BBB = ((tril(A(1:ihi(2),1:ihi(2)),-1)+eye(ihi(2),ihi(2)))'*A(1:ihi(2),ilo(3):ihi(3)));
   BBB = BBB + A(ilo(3):m,1:ihi(2))'*A(ilo(3):m,ilo(3):ihi(3));
   BBB = triu( T(1:ihi(2),1:ihi(2)) )'*BBB;
   A(1:ihi(2),ilo(3):ihi(3)) = A(1:ihi(2),ilo(3):ihi(3)) - (tril(A(1:ihi(2),1:ihi(2)),-1)+eye(ihi(2),ihi(2)))*BBB;
   A(ilo(3):m,ilo(3):ihi(3)) = A(ilo(3):m,ilo(3):ihi(3)) - A(ilo(3):m,1:ihi(2))*BBB;
%
   [ A(ilo(3):m,ilo(3):ihi(3)) ] = geqr2( A(ilo(3):m,ilo(3):ihi(3)) );
%
   [ T(ilo(3):ihi(3),ilo(3):ihi(3)) ] = larft( A(ilo(3):m,ilo(3):ihi(3)) );
%
   T(ilo(1):ihi(2),ilo(3):ihi(3)) = A(ilo(3):ihi(3),ilo(1):ihi(2))' * (tril(A(ilo(3):ihi(3),ilo(3):ihi(3)),-1)+eye(nb(3),nb(3))) ;
   T(ilo(1):ihi(2),ilo(3):ihi(3)) = T(ilo(1):ihi(2),ilo(3):ihi(3)) + A(ihi(3)+1:m,ilo(1):ihi(2))' * A(ihi(3)+1:m,ilo(3):ihi(3))  ;
   T(ilo(1):ihi(2),ilo(3):ihi(3)) = - ( T(ilo(1):ihi(2),ilo(1):ihi(2)) * T(ilo(1):ihi(2),ilo(3):ihi(3)) * T(ilo(3):ihi(3),ilo(3):ihi(3)) ) ;
%
   BBB(ilo(1):ihi(2),1:nb(3)) = triu( T(ilo(1):ihi(2),ilo(1):ihi(2)) )*[A(ilo(3):ihi(3),ilo(1):ihi(2))'];
   BBB(ilo(1):ihi(2),1:nb(3)) = BBB(ilo(1):ihi(2),1:nb(3)) + T(ilo(1):ihi(2),ilo(3):ihi(3))*( tril(A(ilo(3):ihi(3),ilo(3):ihi(3)),-1)+eye(nb(3)) )';
   BBB(ilo(3):ihi(3),1:nb(3)) = triu( T(ilo(3):ihi(3),ilo(3):ihi(3)) )*[tril(A(ilo(3):ihi(3),ilo(3):ihi(3)),-1)+eye(nb(3))]';
%
   Q(ilo(1):ihi(2),ilo(3):ihi(3)) = - (tril(A(ilo(1):ihi(2),ilo(1):ihi(2)),-1)+eye(ihi(2),ihi(2))) * BBB(ilo(1):ihi(2),1:nb(3));
%
   Q(ilo(3):ihi(3),ilo(3):ihi(3)) = ( tril( A(ilo(3):ihi(3),ilo(3):ihi(3)), -1 ) + eye(nb(3),nb(3) ) ) * triu( BBB(ilo(3):ihi(3),1:nb(3)) ) ;
   Q(ilo(3):ihi(3),ilo(3):ihi(3)) = eye(nb(3),nb(3)) - Q(ilo(3):ihi(3),ilo(3):ihi(3));
   Q(ilo(3):ihi(3),ilo(3):ihi(3)) = Q(ilo(3):ihi(3),ilo(3):ihi(3)) - A(ilo(3):ihi(3),ilo(1):ihi(2)) * BBB(ilo(1):ihi(2),1:nb(3));
%
   Q(ihi(3)+1:m,ilo(3):ihi(3)) = - A(ihi(3)+1:m,ilo(3):ihi(3)) * triu( BBB(ilo(3):ihi(3),1:nb(3)) );
   Q(ihi(3)+1:m,ilo(3):ihi(3)) = Q(ihi(3)+1:m,ilo(3):ihi(3)) - A(ihi(3)+1:m,ilo(1):ihi(2)) * BBB(ilo(1):ihi(2),1:nb(3)) ;
%
%
%
   R = triu( A(1:n,1:n) );
%
   fprintf('|| A - Q*R || / || A ||   = %6.1e\n',norm(As(1:m,1:n)-Q(1:m,1:n)*R(1:n,1:n),'fro')/nrmA);
   fprintf('|| I - Q''*Q ||            = %6.1e\n',norm(eye(n) - Q(1:m,1:n)'*Q(1:m,1:n),'fro'));
