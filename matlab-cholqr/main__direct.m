%
   clear
%
   m = 20;
   n = 10;
   log10KA = 5;
%
   U = randn(m,n); [U,~]=qr(U,0);
   V = randn(n,n); [V,~]=qr(V,0);
   S = diag( 10.^( linspace( 0, log10KA, n ) ) );
   A = U * S * V';
   clear U S V;
%
   Q = zeros(m,n);
   R = zeros(n,n);
%
   R = A'*A;
   R = chol( R, 'upper' );
   Q = A / R;
%
   R1 = Q'*Q;
   R1 = chol( R1, 'upper' );
   Q = Q / R1;
   R = R1 * R;
%
   orth = norm(eye(n) - Q'*Q, 'fro');
   repres = norm( A - Q * R, 'fro') / norm( A, 'fro' ) ;
   fprintf('|| I - Q''* Q ||           = %6.1e\n', orth );
   fprintf('|| A - Q * R || / || A || = %6.1e\n', repres );
%
