%
   clear
   clear global
%
   fprintf('\n');
%
   m = 300;
   n = 117;
   log10KA = 11;
%
   mt = 410;
%  nb_lvl = [128, 64, 32, 8, 2 ];
   nb_lvl = [128, 64, 32 ];
%
   n_lvl = size( nb_lvl, 2);
%
   if ( m < n ) fprintf('m < n\n'); return; end
   U = randn(m,n); [U,~]=qr(U,0);
   V = randn(n,n); [V,~]=qr(V,0);
   S = diag( 10.^( linspace( 0, -log10KA, n ) ) );
   A = U * S * V';
   clear U S V
   As = A;
%
   Q = zeros(m,n);
   T = zeros(mt,n);
%  [ A, T, Q ] = lila_geqrf_v05_w00_level1( m, n, 1, mt, A, T, Q );
%  [ A, T, Q ] = lila_geqrf_v05_w01_level1( m, n, 1, mt, A, T, Q );
%  [ A, T, Q ] = lila_geqrf_v05_w02_level1( m, n, 1, mt, A, T, Q );
  [ A, T, Q ] = lila_geqrf_v05_w03_level1( n_lvl, 1, nb_lvl, m, n, 1, mt, A, T, Q );
%  [ A, T, Q ] = lila_geqrf_v05_w04_level1( m, n, 1, mt, A, T, Q );
%
%   [ A, T ] = lila_geqrf_u05_w03_level1( m, n, 1, mt, A, T );
%   [ Q ] = lila_orgqrf_v05_w03( m, n, 1, mt, A, T, Q );
%
   TT = lapack_larft( A );
   V = tril(A(1:m,1:n),-1)+eye(m,n);
   H = eye(m,m) - V * TT * V';
   fprintf('\n||H''*H-I|| = %d',norm(H'*H-eye(m),'fro'))
   fprintf('                 ||tril(H''*As)|| = %d',norm(tril(H'*As,-1),'fro')/norm(A,'fro'))
   fprintf('              ||triu(H''*As)-R|| = %d\n\n',norm(triu(H'*As)-triu(A),'fro')/norm(A,'fro'))
%  fprintf('||qr(As)-A|| = %d',norm(qr(As,0)-A,'fro')/norm(As,'fro'))
   R = triu(A(1:n,1:n));
   fprintf('               ||Q''Q - I|| = %e', norm(Q(1:m,1:n)'*Q(1:m,1:n) - eye(n), 'fro'));
   fprintf('                  ||A - Q*R|| / || A || = %e\n\n', norm(As(1:m,1:n) - Q(1:m,1:n)*R(1:n,1:n), 'fro') / norm(As(1:m,1:n), 'fro'));
%
%   [1:n;T]
