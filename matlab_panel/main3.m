%
   clear
   clear global
%
   m = 23;
   ni = [ 4, 4, 4, 4 ];
   n = sum(ni);
   log10KA = 1;
%
   if ( m < n ) fprintf('m < n\n'); return; end
   U = randn(m,n); [U,~]=qr(U,0);
   V = randn(n,n); [V,~]=qr(V,0);
   S = diag( 10.^( linspace( 0, -log10KA, n ) ) );
   A = U * S * V';
   clear U S V
   As = A;
%
   mt = n;
%
   Q = zeros(m,n);
   T = zeros(n,n);
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

%      [ A, T, Q, D ] = lila_geqr2_w0b_panel( m, n, 1, A, T, Q, D );
%
      [ A, T, Q, D ] = lila_orghr_w0b_panel( m, ni(1), 1, 1, A, T, Q, D );
      [ A, T, Q, D ] = lila_orghr_w0b_panel( m, ni(2), 1, ni(1)+1, A, T, Q, D );
      [ A, T, Q, D ] = lila_orghr_w0b_panel( m, ni(3), ni(1)+ni(2)+1, ni(1)+ni(2)+1, A, T, Q, D );
      [ A, T, Q, D ] = lila_orghr_w0b_panel( m, ni(4), ni(1)+ni(2)+1, ni(1)+ni(2)+ni(3)+1, A, T, Q, D );

	fprintf('n = %d, i = %d, j = %d\n',ni(1), 1, 1);
	fprintf('n = %d, i = %d, j = %d\n',ni(2), 1, ni(1)+1);
	fprintf('n = %d, i = %d, j = %d\n',ni(3), ni(1)+ni(2)+1, ni(1)+ni(2)+1 );
	fprintf('n = %d, i = %d, j = %d\n',ni(4), ni(1)+ni(2)+1, ni(1)+ni(2)+ni(3)+1);
%
%  Checks
%
   R = triu(A(1:n,1:n));
   V = tril(A,-1)+eye(m,n);
%
   T1 = T(1:ni(1)+ni(2),1:ni(1)+ni(2));
   T2 = T(ni(1)+ni(2)+1:n,ni(1)+ni(2)+1:n);
%
   H1 = ( eye(m) - V(1:m,1:ni(1)+ni(2)) * ( T1 * V(1:m,1:ni(1)+ni(2))' ) );
   H2 = ( eye(m) - V(1:m,ni(1)+ni(2)+1:n) * ( T2 * V(1:m,ni(1)+ni(2)+1:n)' ) );
   H  = H1*H2;
%
   fprintf('\n');
   fprintf('|| Q''*Q - I ||            = % 5.3e\n',norm(Q'*Q-eye(n),'fro'))
   fprintf('|| A - Q*R || / || A ||   = % 5.3e\n', norm(As(1:m,1:n) - Q(1:m,1:n)*R(1:n,1:n), 'fro') / norm(As(1:m,1:n), 'fro'));
%
   fprintf('|| H''*H - I || (1,1)      = % 5.3e\n', norm( H(1:m,1:n)'*H(1:m,1:n) - eye(n), 'fro'));
   fprintf('|| H''*H - I || (1,2)      = % 5.3e\n', norm( H(1:m,n+1:m)'*H(1:m,n+1:m) - eye(m-n), 'fro'));
   fprintf('|| H''*H - I || (2,2)      = % 5.3e\n', norm( H(1:m,1:n)'*H(1:m,n+1:m), 'fro'));
%
   fprintf('|| H  -  Q ||             = % 5.3e\n', norm( H(1:m,1:n)-Q, 'fro') );
%
   fprintf('|| H*R - As || / || A ||  = % 5.3e\n', norm( H(1:m,1:n)*R-As, 'fro') / norm(As(1:m,1:n), 'fro'));
   fprintf('|| H''*As - R || / || A || = % 5.3e\n', norm( H(1:m,1:n)'*As-R, 'fro') / norm(As(1:m,1:n), 'fro'));
   fprintf('|| H''*As || / || A ||     = % 5.3e\n', norm( H(1:m,n+1:m)'*As, 'fro') / norm(As(1:m,1:n), 'fro'));
%
   fprintf('\n');
%
