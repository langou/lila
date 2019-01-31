%
   clear
   clear global
%
   fprintf('\n');
%
   m = 17;
   i = 2;
   j = i;
   n = 12;
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
   TT = zeros(n+i-1,n+i-1);
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

   AA = A; QQ = Q;
   [ AA, TT, QQ, DD ] = lila_geqr2_w0b_panel( m, n, i, j, AA, TT, QQ, D );
%
   n1 = ceil((n+i-1)/2);
   n2 = n - n1;
%
   [ A, T, Q, D ] = lila_geqr2_w0b_panel( m, n1, i, j, A, T, Q, D );
   [ A, T, Q, D ] = lila_geqr2_w0b_panel( m, n2, i, n1+j, A, T, Q, D );
%
%  Checks
%
   RR = triu(qr(As(i:m,i:n+i-1),0));
   VV = tril(qr(As(i:m,i:n+i-1),0),-1) + eye(m-i+1,n);
   HH = ( eye(m-i+1) - VV * ( TT(i:n+i-1,i:n+i-1) * VV' ) );
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
