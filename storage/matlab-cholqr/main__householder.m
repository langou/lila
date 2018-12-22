%
   clear
%
   m = 10;
   n = 5;
   log10KA = 3;
%
   U = randn(m,n); [U,~]=qr(U,0);
   V = randn(n,n); [V,~]=qr(V,0);
   S = diag( 10.^( linspace( 0, log10KA, n ) ) );
   A = U * S * V';
   clear U S V;
%

   C = A;
   Asave = A;

   for j = 1:min( [ m n ]),

   T = larft( A(1:m,1:j-1) );
%  [ A(1:m,j) ] = larfb( A(1:m,1:j-1), T, A(1:m,j) );

   A(1:m,j) = A(1:m,j) - (tril(A(1:m,1:j-1),-1)+eye(m,j-1)) * ( T' * ( (tril(A(1:m,1:j-1),-1)+eye(m,j-1))' * A(1:m,j) ) ) ;



%  for i = 1:j-1,
%     [ A(i:m,j) ] = larfL( A(i:m,i), A(i:m,j) ); 
%  end

   [ A(j:m,j) ] = larfg( A(j:m,j) ); 

   end

%  check


   Q = eye(m,m);
   for j = min( [ m n ]):-1:1,
      [ Q(j:m,j:m) ] = larfL( A(j:m,j), Q(j:m,j:m) ); 
   end

  R = triu( A );
  fprintf('|| A - Q*R || / || A || = %6.1e\n',norm(Asave-Q*R)/norm(Asave));
  fprintf('|| I - Q''*Q || = %6.1e\n',norm(eye(m) - Q'*Q));


