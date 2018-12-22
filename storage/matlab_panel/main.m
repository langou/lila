%
   clear
   clear global
%
   fprintf('\n');
%
   m = 15;
   n = 11;
   log10KA = 2;
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
   T = zeros(n,n);
%
%   [ A, T, Q ] = lila_geqr2_w0a( m, n, 1, A, T, Q );
%
%   [ A, T, Q ] = lila_geqr2_w0b( m, n, A, T, Q );  % This script only does a full block with LU 
%
%
%
      Q(1:m,1:n) = A(1:m,1:n);
      D = -zeros(n,1);
      A = rand(m,n);

      [ Q, R ] = Cholesky_qr( Q );
      for ii=1:n,
         for jj=ii:n,
           A(ii,jj) = R(ii,jj);
         end
      end

      AA = A; TT = T; QQ = Q;
      [ AA, TT, QQ, DD ] = lila_geqr2_w0b_panel( m, n, 1, AA, TT, QQ, D );
%      [ A, T, Q, D ] = lila_geqr2_w0b_panel( m, n, 1, A, T, Q, D );
%
      n1 = ceil(n/2);
      n2 = n - n1;
%
      [ A, T, Q, D ] = lila_geqr2_w0b_panel( m, n1, 1, A, T, Q, D );
      [ A, T, Q, D ] = lila_geqr2_w0b_panel( m, n2, n1+1, A, T, Q, D );
%
%
%  Checks
%
   RR = triu(qr(As,0));
   VV = tril(qr(As,0),-1) + eye(m,n);
   TT = lapack_larft( VV );
%
   R = triu(A(1:n,1:n));
   V = tril(A,-1)+eye(m,n);
%
   H = ( eye(m) - V * ( T * V' ) );
   HH = ( eye(m) - V * ( TT * V' ) );
%
   fprintf('\n|| Q''*Q - I || = % 5.3e\n',norm(Q'*Q-eye(n),'fro'))
%
   fprintf('|| H''*H - I || = % 5.3e\n', norm( H(1:m,1:n)'*H(1:m,1:n) - eye(n), 'fro'));
   fprintf('|| H''*H - I || = % 5.3e\n', norm( H(1:m,n+1:m)'*H(1:m,n+1:m) - eye(m-n), 'fro'));
   fprintf('|| H''*H - I || = % 5.3e\n', norm( H(1:m,1:n)'*H(1:m,n+1:m), 'fro'));
%
   fprintf('|| H  -  Q ||  = % 5.3e\n', norm( H(1:m,1:n)-Q, 'fro') );
%
   fprintf('|| H''*As - R || / || A || = % 5.3e\n', norm( H(1:m,1:n)'*As-R, 'fro') / norm(As(1:m,1:n), 'fro'));
   fprintf('|| H''*As || / || A ||     = % 5.3e\n', norm( H(1:m,n+1:m)'*As, 'fro') / norm(As(1:m,1:n), 'fro'));
%
   fprintf('|| A - Q*R || / || A ||   = % 5.3e\n', norm(As(1:m,1:n) - Q(1:m,1:n)*R(1:n,1:n), 'fro') / norm(As(1:m,1:n), 'fro'));
