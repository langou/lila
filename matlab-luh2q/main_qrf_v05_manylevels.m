%
   clear
   clear global
%
   global nb_lvl1;
   global ii_lvl2;
   global nb_lvl2;
%
   fprintf('\n');
   m = 100;
%   mt = 9;
   log10KA = 2;
   nb_lvl2{1} = [   35, 3, 2, 4 ];
   nb_lvl2{2} = [   3, 1, 4 ];
   nb_lvl2{3} = [   1, 2, 5, 2 ];
%  nb_lvl2{4} = [   1, 1  ];
%  nb_lvl2{5} = [   4, 11  ];
%
   nb_blocks_lvl1 =size(nb_lvl2,2);
   nb_lvl1 = zeros(1,nb_blocks_lvl1);
   nb_blocks_lvl2 = zeros(1,nb_blocks_lvl1);
   for i = 1:nb_blocks_lvl1,
      nb_lvl1(i) = sum( nb_lvl2{i} );
      nb_blocks_lvl2(i) = size(nb_lvl2{i},2);
   end
   n = sum(nb_lvl1);
   if ( m < n ) fprintf('m < n\n'); return; end
%
%  we are checking w03
   mt =  11;
%
   U = randn(m,n); [U,~]=qr(U,0);
   V = randn(n,n); [V,~]=qr(V,0);
   S = diag( 10.^( linspace( 0, -log10KA, n ) ) );
   A = U * S * V';
   clear U S V
   As = A;
%
   Q = randn(m,n);
   T = zeros(mt,n);
%  [ A, T, Q ] = lila_geqrf_v05_w00_level1( m, n, 1, mt, A, T, Q );
%  [ A, T, Q ] = lila_geqrf_v05_w01_level1( m, n, 1, mt, A, T, Q );
%  [ A, T, Q ] = lila_geqrf_v05_w02_level1( m, n, 1, mt, A, T, Q );
   [ A, T, Q ] = lila_geqrf_v05_w03_level1( m, n, 1, mt, A, T, Q );
%
   TT = lapack_larft( A );
   V = tril(A(1:m,1:n),-1)+eye(m,n);
   H = eye(m,m) - V * TT * V';
   fprintf('||H''*H-I|| = %d',norm(H'*H-eye(m),'fro'))
   fprintf('                 ||tril(H''*As)|| = %d',norm(tril(H'*As,-1),'fro')/norm(A,'fro'))
   fprintf('              ||triu(H''*As)-R|| = %d\n\n',norm(triu(H'*As)-triu(A),'fro')/norm(A,'fro'))
   fprintf('||qr(As)-A|| = %d',norm(qr(As,0)-A,'fro')/norm(As,'fro'))
   R = triu(A(1:n,1:n));
   fprintf('               ||Q''Q - I|| = %e', norm(Q(1:m,1:n)'*Q(1:m,1:n) - eye(n), 'fro'));
   fprintf('                  ||A - Q*R|| / || A || = %e\n\n', norm(As(1:m,1:n) - Q(1:m,1:n)*R(1:n,1:n), 'fro') / norm(As(1:m,1:n), 'fro'));
%
%   [1:n;T]
