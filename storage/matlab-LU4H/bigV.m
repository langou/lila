%
   clear
%
   nb = [ 3, 4, 5, 4, 9, 2, 3 11, 3, 6, 9, 10, 14, 8];
   m = 181;
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
   As = A;
   nrmA = norm(A,'fro');
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
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
   fprintf('Local check on block 1\n');
   fprintf('||Q''Q - I|| = %f\n', norm(Q'*Q - eye(ihi(1)), 'fro'));
   fprintf('||A - Q*R|| = %f\n\n', norm(As(1:m,1:ihi(1)) - Q*R, 'fro') / norm(As(1:m,1:ihi(1)), 'fro'));
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
   mm = m;
   for k = 2:nb_block,
%
%        Apply V to second block
      for j = 1:ihi(k-1),
         [ A(j:m,ilo(k):ihi(k)) ] = larfL( A(j:m,j), A(j:m,ilo(k):ihi(k)) );
      end
%
%
%
%       start of GEQRF
      for j = ilo(k):ihi(k),
         [ A(j:m,j) ] = larfg( A(j:m,j) );
         [ A(j:m,j+1:ihi(k)) ] = larfL( A(j:m,j), A(j:m,j+1:ihi(k)) );
      end
%        end of GEQRF
%
%
%
%        start of ORGQR
      mm = mm - nb(k-1);
      nn = nb(k);
      AA = A(ilo(k):m,ilo(k):ihi(k));
%
      QQ = zeros( mm, nn );
      QQ(1:nn,1:nn) = eye( nn, nn );
      for j = nn:-1:1,
         [ QQ(j:mm,j:nn) ] = larfL( AA(j:mm,j), QQ(j:mm,j:nn) );
      end
%        end of ORGQR
%
%
%
%        start ORMQR
      QQ = [ zeros(ihi(k-1),nb(k)) ; QQ];
      for j = ihi(k-1):-1:1,
         [ QQ(j:m,1:nn) ] = larfL( A(j:m,j), QQ(j:m,1:nn) );
      end
%         end ORMQR
%
%
%         Local check on kth block
%
      Q(1:m,ilo(k):ihi(k)) = QQ;
      R = triu(A(1:ihi(k),1:ihi(k)));
      fprintf('Local check on block %d\n', k);
      fprintf('||Q''Q - I|| = %f\n', norm(QQ'*QQ - eye(nb(k)), 'fro'));
      fprintf('||A - Q*R|| = %f\n\n', norm(As(1:m,1:ihi(k)) - Q*R, 'fro') / norm(As(1:m,1:ihi(k)), 'fro'));

   end
%
%
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


return


%
%  Apply V to second block
%
   for j = 1:ihi(1),
      [ A(j:m,ilo(2):ihi(2)) ] = larfL( A(j:m,j), A(j:m,ilo(2):ihi(2)) );
   end
%
%  start of GEQRF
   for j = ilo(2):ihi(2),
      [ A(j:m,j) ] = larfg( A(j:m,j) );
      [ A(j:m,j+1:ihi(2)) ] = larfL( A(j:m,j), A(j:m,j+1:ihi(2)) );
   end
%  end of GEQRF
%
%  start of ORGQR
   R2 = triu(A(ilo(2):ihi(2),ilo(2):ihi(2)));
   bR = triu(A(1:ihi(2),1:ihi(2)));
%
   m2 = m - nb(1);
   n2 = nb(2);
   A2 = A(ilo(2):m,ilo(2):ihi(2));
%
   Q2 = zeros( m2, n2 );
   Q2(1:n2,1:n2) = eye( n2, n2 );
   for j = n2:-1:1,
      [ Q2(j:m2,j:n2) ] = larfL( A2(j:m2,j), Q2(j:m2,j:n2) );
   end
%  end of ORGQR
%
%  start ORMQR
   Q2 = [ zeros(nb(1),nb(2)) ; Q2];
   for j = nb(1):-1:1,
      [ Q2(j:m,1:n2) ] = larfL( A(j:m,j), Q2(j:m,1:n2) );
   end
%  end ORMQR
%
   Q(1:m,ilo(2):ihi(2)) = Q2;
%  Local check on second block
%
   fprintf('Local check on block two\n');
   fprintf('||Q''Q - I|| = %f\n', norm(Q2'*Q2 - eye(nb(2)), 'fro'));
   fprintf('||A - Q*R|| = %f\n\n', norm(As(1:m,1:ihi(2)) - Q*bR, 'fro') / norm(As(1:m,ilo(2):ihi(2)), 'fro'));
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Apply V to third block
%
   for j = 1:ihi(2),
      [ A(j:m,ilo(3):ihi(3)) ] = larfL( A(j:m,j), A(j:m,ilo(3):ihi(3)) );
   end
%
%  start of GEQRF
   for j = ilo(3):ihi(3),
      [ A(j:m,j) ] = larfg( A(j:m,j) );
      [ A(j:m,j+1:ihi(3)) ] = larfL( A(j:m,j), A(j:m,j+1:ihi(3)) );
   end
%  end of GEQRF
%
%  start of ORGQR
   R = triu(A(1:ihi(3),1:ihi(3)));
%
   m3 = m2 - nb(2);
   n3 = nb(3);
   A3 = A(ilo(3):m,ilo(3):ihi(3));
%
   Q3 = zeros( m3, n3 );
   Q3(1:n3,1:n3) = eye( n3, n3 );
   for j = n3:-1:1,
      [ Q3(j:m3,j:n3) ] = larfL( A3(j:m3,j), Q3(j:m3,j:n3) );
   end
%  end of ORGQR
%
%  start ORMQR
   Q3 = [ zeros(ihi(2),nb(3)) ; Q3];
   for j = ihi(2):-1:1,
      [ Q3(j:m,1:n3) ] = larfL( A(j:m,j), Q3(j:m,1:n3) );
   end
%  end ORMQR
%
   Q(1:m,ilo(3):ihi(3)) = Q3;
%  Local check on second block
%
   fprintf('Local check on block three\n');
   fprintf('||Q''Q - I|| = %f\n', norm(Q'*Q - eye(n), 'fro'));
   fprintf('||A - Q*R|| = %f\n\n', norm(As(1:m,1:ihi(3)) - Q*R, 'fro') / norm(As(1:m,ilo(3):ihi(3)), 'fro'));

%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

