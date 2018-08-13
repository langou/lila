%
   clear
%
   m = 20;
   n = 10;
   log10KA = 3;
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
   for j=1:n,

      A1 = A(1:m,1:j-1)'*A(1:m,1:j-1);
      T1 = inv( A1 );
      T2 = 2*T1 - T1*A1*T1;

      w = (eye(m) - A(1:m,1:j-1) * T2 * A(1:m,1:j-1)' ) * A(1:m,j);
      Q(:,j) = w / norm(w);

   end
   R = Q'*A;
%
   orth = norm(eye(n) - Q'*Q, 'fro');
   repres = norm( A - Q * R, 'fro') / norm( A, 'fro' ) ;
   fprintf('|| I - Q''* Q ||           = %6.1e\n', orth );
   fprintf('|| A - Q * R || / || A || = %6.1e\n', repres );
%
