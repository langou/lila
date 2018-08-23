%
   clear
%
   nb = [ 3, 4, 5, 4, 9, 2, 3 11, 3, 6, 9, 10, 14, 8, 7, 12, 14, 9, 2, 5];
   m = 281;
%
   n = sum(nb);
   nb_block = size(nb,2);
%
   log10KA = 10;
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
   As = A;
   nrmA = norm(A,'fro');
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
   V = tril( A(1:m,1:ihi(1)), -1 ) + eye(m,ihi(1));
%
%  ORGQR on the first block
%
   R = triu(A(1:ihi(1),1:ihi(1)));
   Q = eye( m, ihi(1) );
   for j = ihi(1):-1:1,
      [ Q(j:m,j:ihi(1)) ] = larfL( A(j:m,j), Q(j:m,j:ihi(1)) );
   end
%
%  Local check on first block
%
   fprintf('Local check on block 1');
   fprintf('               Semi-local check on block 1\n');
   fprintf('||Q''Q - I|| = %f', norm(Q'*Q - eye(nb(1)), 'fro'));
   fprintf('                ||A - Q*R|| = %f\n\n', norm(As(1:m,1:ihi(1)) - Q*R, 'fro') / norm(As(1:m,1:ihi(1)), 'fro'));
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
   ml= m;
   for k = 2:nb_block,
%
      ml = ml - nb(k-1);
      nl = nb(k);
%     ORMQRf
%
%     for j = 1:ihi(k-1),
%        [ A(j:m,ilo(k):ihi(k)) ] = larfL( A(j:m,j), A(j:m,ilo(k):ihi(k)) );
%     end
%
      [ A ] = lila_ormqrf_v0( m, ihi(k-1), nb(k), A, 1, 1, lda, A, 1, ilo(k), lda );
%
%
%
%     start of GEQRF
%     for j = ilo(k):ihi(k),
%        [ A(j:m,j) ] = larfg( A(j:m,j) );
%        [ A(j:m,j+1:ihi(k)) ] = larfL( A(j:m,j), A(j:m,j+1:ihi(k)) );
%     end
%
      [ A ] = lila_geqrf_v0( ml, nl, A, ilo(k), ilo(k), lda );
%
%     end of GEQRF
%
%
%
%     start of ORGQR
      nn = nb(k);
      AA = A(ilo(k):m,ilo(k):ihi(k));
%
      QQ = zeros( ml, nn );
      QQ(1:nn,1:nn) = eye( nn, nn );
      for j = nn:-1:1,
         [ QQ(j:ml,j:nn) ] = larfL( AA(j:ml,j), QQ(j:ml,j:nn) );
      end
%
%    Q(ilo(k):m,ilo(k):ihi(k)) = A(ilo(k):m,ilo(k):ihi(k));
%    [ Q ] = lila_orgqr_v0( ml, nl, Q, ilo(k), ilo(k), lda );
%
%     end of ORGQR
%
%
%
%     start ORMQRbz
      QQ = [ zeros(ihi(k-1),nb(k)) ; QQ];
      for j = ihi(k-1):-1:1,
         [ QQ(j:m,1:nn) ] = larfL( A(j:m,j), QQ(j:m,1:nn) );
      end
%     end ORMQRbz
%
%     [ A ] = lila_ormqrbz_v0( m, ihi(k-1), nb(k), A, 1, 1, lda, Q, 1, ilo(k), ldq );
%
%
%
%         Local check on kth block
%
      Q(1:m,ilo(k):ihi(k)) = QQ;
      R = triu(A(1:ihi(k),1:ihi(k)));
      fprintf('Local check on block %d', k);
      fprintf('               Semi-local check on blocks 1-%d\n', k);
      fprintf('||Q''Q - I|| = %f', norm(QQ'*QQ - eye(nb(k)), 'fro'));
      fprintf('                ||A - Q*R|| = %f\n\n', norm(As(1:m,1:ihi(k)) - Q*R, 'fro') / norm(As(1:m,1:ihi(k)), 'fro'));

   end
%
