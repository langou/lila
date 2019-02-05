%
   clear
   clear global
%
   fprintf('\n');
%
   m = 847;
   i = 158;
   j = i;
   n = 473;
   mt = 34;
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
   T = zeros(mt,n+i-1);
%
   Q(i:m,i:n+i-1) = A(i:m,i:n+i-1);
   D = -zeros(n+i-1,1);
   A = rand(m,n+i-1);
 
   [ Q(i:m,i:n+i-1), R ] = Cholesky_qr( Q(i:m,i:n+i-1) );
   for ii=i:n+i-1,
      for jj=ii:n+i-1,
         A(ii,jj) = R(ii-i+1,jj-i+1);
      end
   end
%
   vb = mt - mod(i-1,mt);
   if (vb > n+i-1), vb = n+i-1; end
%
   while( vb~=0 ),

      [ A, T, Q, D ]     = lila_geqr2_w0b_panel_mt( m, vb, i, j, mt, A, T, Q, D );
      j = j+vb;
      if( j+mt-1 >= n+i-1) vb = n + i -1 + 1 - j; else, vb = mt; end

   end

%
%  Checks
%
   R = triu(A(i:n+i-1,i:n+i-1));
   V = tril(A(i:m,i:n+i-1),-1)+eye(m-i+1,n);
   H = eye(m-i+1);
%
   j = i;
   k = 1;
   vb = mt - mod(i-1,mt);
   if (vb > n+i-1), vb = n+i-1; end
   while( vb ~= 0 ),

      it = mod(j-1,mt)+1;
      H = H * ( eye(m-i+1) - V(1:m-i+1,k:k+vb-1) * ( T(it:it+vb-1,j:j+vb-1) * V(1:m-i+1,k:k+vb-1)' ) ) ;
      j = j+vb;
      k = k+vb;
      if( j+mt-1 >= n+i-1) vb = n + i -1 + 1 - j; else, vb = mt; end

   end

   fprintf('\n');
   fprintf('|| Q''*Q - I ||            = % 5.3e\n',norm(Q(i:m,i:n+i-1)'*Q(i:m,i:n+i-1)-eye(n),'fro'))
   fprintf('|| A - Q*R || / || A ||   = % 5.3e\n', norm(As(i:m,i:n+i-1) - Q(i:m,i:n+i-1)*R(1:n,1:n), 'fro') / norm(As(i:m,i:n+i-1), 'fro'));
%
   fprintf('|| H''*H - I || (1,1)      = % 5.3e\n', norm( H(1:m-i+1,1:n)'*H(1:m-i+1,1:n) - eye(n), 'fro'));
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
