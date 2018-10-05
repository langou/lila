%
   clear
   clear global
%
   global nb_lvl1;
   global ii_lvl2;
   global nb_lvl2;
%
   nb_blocks_lvl1 = 3;
   log10KA = 2;
   nb_lvl2 = cell(1,nb_blocks_lvl1);
%  for i = 1:nb_blocks_lvl1, 
%     nb_blocks_lvl2 = 1;
%     for j = 1:nb_blocks_lvl2, 
%        nb_lvl2{i}(1,j) = 6;
%     end
%  end
   nb_lvl2{1} = [ 4 3 11 4 ];
   nb_lvl2{2} = [ 9 5 7 5 8 ];
   nb_lvl2{3} = [ 6 10 3 4 2 ];
   mt = 10;
%
   nb_lvl1 = zeros(1,nb_blocks_lvl1);
   nb_blocks_lvl2 = zeros(1,nb_blocks_lvl1);
   for i = 1:nb_blocks_lvl1,
      nb_lvl1(i) = sum( nb_lvl2{i} );
      nb_blocks_lvl2(i) = size(nb_lvl2{i},2);
   end
   n = sum(nb_lvl1);
%
   m = 85;
%
   U = randn(m,n); [U,~]=qr(U,0);
   V = randn(n,n); [V,~]=qr(V,0);
   S = diag( 10.^( linspace( 0, -log10KA, n ) ) );
   A = U * S * V';
   clear U S V
   As = A;
%
   Q = zeros(m,n);
   T = zeros(mt,n);
%  [ A, T, Q ] = lila_geqrf_manylevels_level1_w00( m, n, 1, mt, A, T, Q );
%   [ A, T, Q ] = lila_geqrf_manylevels_level1_w01_flopsave( m, n, 1, mt, A, T, Q );
%  [ A, T, Q ] = lila_geqrf_manylevels_level1_w02_flopsave( m, n, 1, mt, A, T, Q );
  [ A, T, Q ] = lila_geqrf_manylevels_level1_w03_flopsave( m, n, 1, mt, A, T, Q );
%
%
   TT = lapack_larft( A );
   V = tril(A(1:m,1:n),-1)+eye(m,n);
   H = eye(m,m) - V * TT * V';
   R = triu(A(1:n,1:n));
%
   fprintf('\n');
   fprintf('||H''H - I|| = %e', norm(H'*H-eye(m),'fro'));
   fprintf('    ||tril(H''As)|| = %e', norm(tril(H'*As,-1),'fro')/norm(A,'fro'));
   fprintf('    ||triu(H''As) - R|| = %e', norm(triu(H'*As)-triu(A),'fro')/norm(A,'fro'));
   fprintf('    ||Q''Q - I|| = %e', norm(Q(1:m,1:n)'*Q(1:m,1:n) - eye(n), 'fro'));
   fprintf('    ||A - QR|| = %e', norm(As(1:m,1:n) - Q(1:m,1:n)*R(1:n,1:n), 'fro') / norm(As(1:m,1:n), 'fro'));
   fprintf('\n');

