%
   clear
%
   m = 10;
   n = 10;
   log10KA = 1.5;
   i = 3;
   nl = 5;
   mt = 10;
%
   if (nl+i-1>n), fprintf('error nl+i-1>n\n'); end
   if ( m < n ), fprintf('m < n\n'); return; end
%
   U = randn(m,n); [U,~]=qr(U,0);
   V = randn(n,n); [V,~]=qr(V,0);
   S = diag( 10.^( linspace( 0, -log10KA, n ) ) );
   A = U * S * V';
   clear U S V
   As = A;
   nrmA = norm(As,'fro');
%
   Q = randn(m,n);
   T = randn(mt,n);
%
   [ A, T, Q ] = lila_geqr2_v05_q00(m, nl, i, mt, A, T, Q);
%
   ilo = i;
   ihi = i+nl-1;
   ml = m-i+1;
%
   R = triu(A(ilo:ihi,ilo:ihi));
   TT = lapack_larft( A(i:m,ilo:ihi) );
   V = tril(A(i:m,ilo:ihi),-1)+eye(ml,nl);
   H = eye(ml,ml) - V * TT * V';
%  fprintf('||qr(As)-A|| = %e',norm(qr(As(i:m,ilo:ihi),0)-A(i:m,ilo:ihi),'fro')/norm(As,'fro'))
   fprintf('       ||Q''Q - I|| = %e', norm(Q(i:m,ilo:ihi)'*Q(i:m,ilo:ihi) - eye(nl), 'fro'));
   fprintf('       ||A - Q*R|| = %e', norm(As(i:m,ilo:ihi) - Q(i:m,ilo:ihi)*R, 'fro'));
   fprintf('       ||H''*H-I|| = %e',norm(H'*H-eye(ml),'fro'))
   fprintf('       ||tril(H''*As)|| = %e',norm(tril(H'*As(i:m,ilo:ihi),-1),'fro')/nrmA);
   fprintf('       ||triu(H''*As)-R|| = %e\n',norm(triu(H'*As(i:m,ilo:ihi))-triu(A(i:m,ilo:ihi)),'fro')/nrmA);
