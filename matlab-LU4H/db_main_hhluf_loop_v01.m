%
   clear
%
   nb = [ 3, 4, 5, 4, 9, 2, 3 11, 3, 6, 9, 10, 14, 8, 7, 12, 14, 9, 2, 5];
   m = 281;
%
   n = sum(nb);
   nb_block = size(nb,2);
%
   log10KA = 2;
%
   U = randn(m,n); [U,~]=qr(U,0);
   V = randn(n,n); [V,~]=qr(V,0);
   S = diag( 10.^( linspace( 0, -log10KA, n ) ) );
   A = U * S * V';
   clear U S V;
%
   ilo = zeros(nb_block,1);
   ilo(1) = 1;
   for i = 2:nb_block,
      ilo(i) = ilo(i-1) + nb(i-1);
   end
%
   ihi = zeros(nb_block,1);
   ihi = nb(1);
   for i = 2:nb_block,
      ihi(i) = ihi(i-1) + nb(i);
   end
%
   T = zeros(n,n);
   As = A;
   nrmA = norm(A,'fro');
   Q = zeros( m, n );
%
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%  GEQRF on the first block
%
   for j = 1:ihi(1), 
      [ A(j:m,j) ] = larfg( A(j:m,j) );
      [ A(j:m,j+1:ihi(1)) ] = larfL( A(j:m,j), A(j:m,j+1:ihi(1)) );
   end
%
   T(1:ihi(1),1:ihi(1)) = larft( A(1:m,1:ihi(1)) );
%
%  ORGQR on the first block
%
   R = triu(A(1:ihi(1),1:ihi(1)));
   Q(1:m,ilo(1):ihi(1)) = eye( m, ihi(1) );
   for j = ihi(1):-1:1,
      [ Q(j:m,j:ihi(1)) ] = larfL( A(j:m,j), Q(j:m,j:ihi(1)) );
   end
%
%  Local check on first block
%
   fprintf('Local check on block 1');
   fprintf('               Semi-local check on block 1\n');
   fprintf('||Q''Q - I|| = %f', norm(Q(1:m,ilo(1):ihi(1))'*Q(1:m,ilo(1):ihi(1)) - eye(nb(1)), 'fro'));
   fprintf('                ||A - Q*R|| = %f\n\n', norm(As(1:m,1:ihi(1)) - Q(1:m,ilo(1):ihi(1))*R(ilo(1):ihi(1),ilo(1):ihi(1)), 'fro') / norm(As(1:m,1:ihi(1)), 'fro'));
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%
%                       Looooooooopping time
%
   lda = -1;
   ldq = -1;
   ml= m;
   for k = 2:nb_block-1,
%
      ml = ml - nb(k-1);
      nl = nb(k);
%
%     start of ORMQRF
%     for j = 1:ihi(k-1),
%        [ A(j:m,ilo(k):ihi(k)) ] = larfL( A(j:m,j), A(j:m,ilo(k):ihi(k)) );
%     end
%     start of ORMQRF
%
%     [ A ] = lila_ormqrf_v00( m, ihi(k-1), nb(k), A, 1, 1, lda, A, 1, ilo(k), lda );
     [ A ] = lila_ormqrf_v01( m, ihi(k-1), nb(k), A, 1, 1, lda, A, 1, ilo(k), lda, T );
%
%     start of GEQRF
%     for j = ilo(k):ihi(k),
%        [ A(j:m,j) ] = larfg( A(j:m,j) );
%        [ A(j:m,j+1:ihi(k)) ] = larfL( A(j:m,j), A(j:m,j+1:ihi(k)) );
%     end
%     end of GEQRF
%
      [ A ] = lila_geqrf_v00( ml, nl, A, ilo(k), ilo(k), lda );
%
%
      T(ilo(k):ihi(k),ilo(k):ihi(k) ) = larft( A(ilo(k):m,ilo(k):ihi(k) ) );
      T(1:ihi(k-1),ilo(k):ihi(k)) = -( T(1:ihi(k-1),1:ihi(k-1))*A(ilo(k):m,1:ihi(k-1))' )*( (tril(A(ilo(k):m,ilo(k):ihi(k)),-1)+eye(m-ilo(k)+1,nb(k)))*T(ilo(k):ihi(k),ilo(k):ihi(k) ) );
%
%
%     start of ORGQR
%     nn = nb(k);
%     AA = A(ilo(k):m,ilo(k):ihi(k));
%     QQ = zeros( ml, nn );
%     QQ(1:nn,1:nn) = eye( nn, nn );
%     for j = nn:-1:1,
%        [ QQ(j:ml,j:nn) ] = larfL( AA(j:ml,j), QQ(j:ml,j:nn) );
%     end
%     Q(ilo(k):m,ilo(k):ihi(k)) = QQ;
%     end of ORGQR
%
      Q(ilo(k):m,ilo(k):ihi(k)) = A(ilo(k):m,ilo(k):ihi(k));
      [ Q ] = lila_orgqr_v00( ml, nl, Q, ilo(k), ilo(k), ldq );
%
%     start ORMQRbz
%     QQ = [ zeros(ihi(k-1),nb(k)) ; QQ];
%     for j = ihi(k-1):-1:1,
%        [ QQ(j:m,1:nn) ] = larfL( A(j:m,j), QQ(j:m,1:nn) );
%     end
%     end ORMQRbz
%
      [ Q ] = lila_ormqrbz_v00( m, ihi(k-1), nb(k), A, 1, 1, lda, Q, 1, ilo(k), ldq );
%
%
      R = triu(A(1:ihi(k),1:ihi(k)));
      fprintf('||Q''Q - I|| = %e', norm(Q(1:m,1:ihi(k))'*Q(1:m,1:ihi(k)) - eye(ihi(k)), 'fro'));
      fprintf('                ||A - Q*R|| = %e\n\n', norm(As(1:m,1:ihi(k)) - Q(1:m,1:ihi(k))*R(1:ihi(k),1:ihi(k)), 'fro') / norm(As(1:m,1:ihi(k)), 'fro'));

   end
%%%
%
%  Taking out the final block to work through the code on
%
%
      k = nb_block;
      ml = ml - nb(nb_block-1);
      nl = nb(nb_block);
%
      [ A ] = lila_ormqrf_v01( m, ihi(nb_block-1), nb(nb_block), A, 1, 1, lda, A, 1, ilo(nb_block), lda, T );
%
%     Not doing anything with geqrf yet
%
      [ A ] = lila_geqrf_v00( ml, nl, A, ilo(nb_block), ilo(nb_block), lda );
%
      T(ilo(nb_block):ihi(nb_block),ilo(nb_block):ihi(nb_block) ) = larft( A(ilo(nb_block):m,ilo(nb_block):ihi(nb_block) ) );
      T(1:ihi(nb_block-1),ilo(nb_block):ihi(nb_block)) = -( T(1:ihi(nb_block-1),1:ihi(nb_block-1))*A(ilo(nb_block):m,1:ihi(nb_block-1))' )*( (tril(A(ilo(nb_block):m,ilo(nb_block):ihi(nb_block)),-1)+eye(m-ilo(nb_block)+1,nb(nb_block)))*T(ilo(nb_block):ihi(nb_block),ilo(nb_block):ihi(nb_block) ) );
%
%
      Q(ilo(nb_block):m,ilo(nb_block):ihi(nb_block)) = A(ilo(nb_block):m,ilo(nb_block):ihi(nb_block));
      [ Q ] = lila_orgqr_v00( ml, nl, Q, ilo(k), ilo(k), ldq );
%
      [ Q ] = lila_ormqrbz_v00( m, ihi(nb_block-1), nb(nb_block), A, 1, 1, lda, Q, 1, ilo(nb_block), ldq );
%
%
      R = triu(A(1:ihi(nb_block),1:ihi(nb_block)));
      fprintf('||Q''Q - I|| = %e', norm(Q(1:m,1:ihi(nb_block))'*Q(1:m,1:ihi(nb_block)) - eye(ihi(nb_block)), 'fro'));
      fprintf('                ||A - Q*R|| = %e\n\n', norm(As(1:m,1:ihi(nb_block)) - Q(1:m,1:ihi(nb_block))*R(1:ihi(nb_block),1:ihi(nb_block)), 'fro') / norm(As(1:m,1:ihi(nb_block)), 'fro'));
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Checks for the good T
%
   V = ( tril( A(1:m,1:ihi(nb_block)), -1 ) + eye(m,ihi(nb_block)) ) ;
   T1 = larft( V );
   H = (eye(m,m) - V * ( T(1:ihi(nb_block),1:ihi(nb_block)) * ( V' ) ) );
   H1 = (eye(m,m) - V * ( T1 * ( V' ) ) );

   fprintf('||H''*H - I|| = %e', norm(H'*H - eye(m,m), 'fro') );
   fprintf('   ||H''*Q - I|| = %e', norm(H'*Q - eye(m,n), 'fro') );
   fprintf('   ||triu(H''*A) / norm(R)|| = %e\n', norm( H'*As - [ R; zeros(m-n,n) ], 'fro') / norm( As, 'fro') );
   fprintf('\n||H1 - H|| = %e\n', norm(H1 - H, 'fro') );
