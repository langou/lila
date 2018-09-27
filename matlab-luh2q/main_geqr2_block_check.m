%
   fprintf('\n');
   clear
%
   m = 31;
   n =  21;
   log10KA = 1.5;
   i = 3;
   nl = 13;
   mt = 5;
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
%  T = zeros(mt,n);
%
   [ A, T, Q ] = lila_geqr2_v05_q00(m, nl, i, mt, A, T, Q);
%
   ilo = i;
   ihi = i+nl-1;
   ml = m-i+1;
%
   R = triu(A(ilo:ihi,ilo:ihi));
%  TT = lapack_larft( A(i:m,ilo:ihi) );  V = tril(A(i:m,ilo:ihi),-1)+eye(ml,nl);  HH = eye(ml,ml) - V * TT * V';
   fprintf('(1)');
   fprintf('       ||Q''Q - I|| = %e', norm(Q(i:m,ilo:ihi)'*Q(i:m,ilo:ihi) - eye(nl), 'fro'));
   fprintf('       ||A - Q*R|| = %e', norm(As(i:m,ilo:ihi) - Q(i:m,ilo:ihi)*R, 'fro'));
%  fprintf('       ||H''*H-I|| = %e',norm(HH'*HH-eye(ml),'fro'))
%  fprintf('       ||tril(H''*As)|| = %e',norm(tril(HH'*As(i:m,ilo:ihi),-1),'fro')/nrmA);
%  fprintf('       ||triu(H''*As)-R|| = %e',norm(triu(HH'*As(i:m,ilo:ihi))-triu(A(i:m,ilo:ihi)),'fro')/nrmA);
   fprintf('\n');
%
   Q = randn(m,n);
   [ Q ] = lila_orgqrf_v05_w03( m, nl, i, mt, A, T, Q );
   fprintf('(2)');
   fprintf('       ||Q''Q - I|| = %e', norm(Q(i:m,ilo:ihi)'*Q(i:m,ilo:ihi) - eye(nl), 'fro'));
   fprintf('       ||A - Q*R|| = %e', norm(As(i:m,ilo:ihi) - Q(i:m,ilo:ihi)*R, 'fro'));
%
   Q = randn(m,n);
   H = randn(m,ml);
   H(i:m,1:ml) = eye(ml,ml);
   A = [ A H ];
   [ A ] = lila_ormqrf_v05_w03( m, ml, nl, i, n+1, mt, A, T );
   H = A(i:m,n+1:n+ml);
   A = A(1:m,1:n);
   fprintf('       ||H''*H-I|| = %e',norm(H'*H-eye(ml),'fro'))
   fprintf('       ||tril(H''*As)|| = %e',norm(tril(H*As(i:m,ilo:ihi),-1),'fro')/nrmA);
   fprintf('       ||triu(H''*As)-R|| = %e',norm(triu(H*As(i:m,ilo:ihi))-triu(A(i:m,ilo:ihi)),'fro')/nrmA);
   fprintf('\n');
%
   Q = randn(m,n);
   A = [ A [ rand(i-1,nl); As(i:m,ilo:ihi)] ];
   [ A ] = lila_ormqrf_v05_w03( m, nl, nl, i, n+1, mt, A, T );
   R_ = A(i:m,n+1:n+nl);
   A = A(1:m,1:n);
   fprintf('(3)');
   fprintf('       ||tril(H''*As)|| = %e',norm(tril(R_(1:ml,1:nl),-1)/nrmA));
   fprintf('   ||triu(H''*As)-R|| = %e',norm((triu(R_(1:nl,1:nl))-R)/nrmA));
   fprintf('\n');


