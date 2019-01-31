%
   clear
   clear global
%
   m = 47;
   i = 7;
   n = 26;
   nb = 4;
   mt = n+i-1;
   log10KA = 1;
%
   if ( m < n ) fprintf('m < n\n'); return; end
   U = randn(m,n+i-1); [U,~]=qr(U,0);
   V = randn(n+i-1,n+i-1); [V,~]=qr(V,0);
   S = diag( 10.^( linspace( 0, -log10KA, n+i-1 ) ) );
   A = U * S * V';
   clear U S V
   As = A;
%
   T = zeros(n+i-1,n+i-1);
   TT = zeros(mt,n+i-1);
%
   Q(1:m,1:n+i-1) = A(1:m,1:n+i-1);
   D = -zeros(n+i-1,1);
   A = rand(m,n+i-1);

   [ Q, R ] = Cholesky_qr( Q(i:m,i:n+i-1) );
   for ii=i:n+i-1,
      for jj=ii:n+i-1,
         A(ii,jj) = R(ii-i+1,jj-i+1);
      end
   end

   for ii = 1:i-1,
      Q = [zeros(1,n+ii-1);Q];
      Q = [zeros(m-i+1+ii,1),Q];
   end
   for ii = 1:i-1,
      Q(ii,ii) = 1.0e+00;
   end

   j = i;
   if( nb > n+i-1 ), vb = n+i-1; else vb = nb; end
   [ A, T, Q, D ] = lila_orghr_w0b_panel( m, vb, i, j, A, T, Q, D );
   j = j+vb;
   if( j+nb > n+i-1 ) vb = n+i-1-j; else vb = nb; end
   while( vb~= 0 )
      [ A, T, Q, D ] = lila_orghr_w0b_panel( m, vb, i, j, A, T, Q, D );
      j = j+vb;
      if( j+nb > n+i-1 ) vb = n+i-1-j; else vb = nb; end
   end

%
%  Checks
%
   R = triu(A(i:n+i-1,i:n+i-1));
   V = tril(A(i:m,i:n+i-1),-1)+eye(m-i+1,n);
   H = ( eye(m-i+1) - V * ( T(i:n+i-1,i:n+i-1) * V' ) );

%
   fprintf('\n');
   fprintf('|| Q''*Q - I ||            = % 5.3e\n',norm(Q(i:m,i:n+i-1)'*Q(i:m,i:n+i-1)-eye(n),'fro'))
   fprintf('|| A - Q*R || / || A ||   = % 5.3e\n', norm(As(i:m,i:n+i-1) - Q(i:m,i:n+i-1)*R(1:n,1:n), 'fro') / norm(As(i:m,i:n+i-1), 'fro'));
%
   fprintf('|| H''*H - I || (1,1)      = % 5.3e\n', norm( H(1:m-i+1,1:n)'*H(1:m-i+1,1:n) - eye(n), 'fro'));
%   fprintf('|| H''*H - I || (1,2)      = % 5.3e\n', norm( H(1:m-i+1,n+i:m-i+1)'*H(1:m-i+1,n+i:m-i+1) - eye(m-i+1-n), 'fro'));
   fprintf('||   H''*H   || (2,2)      = % 5.3e\n', norm( H(1:m-i+1,1:n)'*H(1:m-i+1,n+i:m-i+1), 'fro'));
%
   fprintf('|| H  -  Q ||             = % 5.3e\n', norm( H(1:m-i+1,1:n)-Q(i:m,i:n+i-1), 'fro') );
%
   fprintf('|| H*R - As || / || A ||  = % 5.3e\n', norm( H(1:m-i+1,1:n)*R-As(i:m,i:n+i-1), 'fro') / norm(As(i:m,i:n+i-1), 'fro'));
   fprintf('|| H''*As - R || / || A || = % 5.3e\n', norm( H(1:m-i+1,1:n)'*As(i:m,i:n+i-1)-R, 'fro') / norm(As(i:m,i:n+i-1), 'fro'));
   fprintf('|| H''*As || / || A ||     = % 5.3e\n', norm( H(1:m-i+1,n+i:m-i+1)'*As(i:m,i:n+i-1), 'fro') / norm(As(i:m,i:n+i-1), 'fro'));
%
   fprintf('\n');
%
