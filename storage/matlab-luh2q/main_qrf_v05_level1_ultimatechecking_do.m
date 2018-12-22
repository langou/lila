%
function [ check ] = main_qrf_v05_level1_ultimatechecking_do( m, n, mt, log10KA  )
%
   global nb_lvl1;
   global ii_lvl2;
   global nb_lvl2;
%
   nb_blocks_lvl1 =size(nb_lvl2,2);
   nb_lvl1 = zeros(1,nb_blocks_lvl1);
   nb_blocks_lvl2 = zeros(1,nb_blocks_lvl1);
   for i = 1:nb_blocks_lvl1,
      nb_lvl1(i) = sum( nb_lvl2{i} );
      nb_blocks_lvl2(i) = size(nb_lvl2{i},2);
   end

   U = randn(m,n); [U,~]=qr(U,0);
   V = randn(n,n); [V,~]=qr(V,0);
   S = diag( 10.^( linspace( 0, -log10KA, n ) ) );
   A = U * S * V';
   clear U S V
   As = A;
%
   Q = randn(m,n);
   T = zeros(mt,n);
%  [ A, T, Q ] = lila_geqrf_manylevels_level1_w00( m, n, 1, mt, A, T, Q );
%  [ A, T, Q ] = lila_geqrf_manylevels_level1_w01( m, n, 1, mt, A, T, Q );
  [ A, T, Q ] = lila_geqrf_manylevels_level1_w02( m, n, 1, mt, A, T, Q );
%   [ A, T, Q ] = lila_geqrf_manylevels_level1_w03( m, n, 1, mt, A, T, Q );
%
   TT = lapack_larft( A );
   V = tril(A(1:m,1:n),-1)+eye(m,n);
   H = eye(m,m) - V * TT * V';
   check = ones(1,6);
   check(1) = norm(H'*H-eye(m),'fro');
   check(2) = norm(tril(H'*As,-1),'fro')/norm(A,'fro');
   check(3) = norm(triu(H'*As)-triu(A),'fro')/norm(A,'fro');
   check(4) = norm(qr(As,0)-A,'fro')/norm(qr(As),'fro');
   R = triu(A(1:n,1:n));
   check(5) = norm(Q(1:m,1:n)'*Q(1:m,1:n) - eye(n), 'fro');
   check(6) = norm(As(1:m,1:n) - Q(1:m,1:n)*R(1:n,1:n), 'fro') / norm(As(1:m,1:n), 'fro');


%triu(H'*As)-triu(A)
%

end
