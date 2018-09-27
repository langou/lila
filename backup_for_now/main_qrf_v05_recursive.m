%
   clear
%
   fprintf('\n');
   m = 100;
   n = 30;
   mt = 10;
   log10KA = 2;
%
   U = randn(m,n); [U,~]=qr(U,0);
   V = randn(n,n); [V,~]=qr(V,0);
   S = diag( 10.^( linspace( 0, -log10KA, n ) ) );
   A = U * S * V';
   clear U S V;
   As = A;
%
   Q = randn(m,n);
   T = zeros(mt,n);
   [ A, T, Q ] = lila_geqrf_v05_recursive( m, n, 1, mt, A, T, Q );
%
   TT = lapack_larft( A );
   V = tril(A(1:m,1:n),-1)+eye(m,n);
   H = eye(m,m) - V * TT * V';
   R = triu(A(1:n,1:n));
   fprintf('||H''*H-I|| = %d',norm(H'*H-eye(m),'fro'))
   fprintf('                 ||tril(H''*As)|| = %d',norm(tril(H'*As,-1),'fro')/norm(A,'fro'))
   fprintf('              ||triu(H''*As)-R|| = %d\n\n',norm(triu(H'*As)-triu(A),'fro')/norm(A,'fro'))
   fprintf('||qr(As)-A|| = %d',norm(qr(As,0)-A,'fro')/norm(As,'fro'))
   fprintf('               ||Q''Q - I|| = %e', norm(Q(1:m,1:n)'*Q(1:m,1:n) - eye(n), 'fro'));
   fprintf('                  ||A - Q*R|| / || A || = %e\n\n', norm(As(1:m,1:n) - Q(1:m,1:n)*R(1:n,1:n), 'fro') / norm(As(1:m,1:n), 'fro'));
%
