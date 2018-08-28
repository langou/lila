%
   clear
%
%   nb = [ 3, 4, 5, 4, 9, 2, 3 11, 3, 6, 9, 10, 14, 8, 7, 12, 14, 9, 2, 5];
   nb = [ 13, 4, 9, 2, 7];
   ib = [ 4, 7, 5, 3, 4, 6, 3, 2, 1 ];
   m = 81;
%
   n = sum(nb);
   nb_block = size(nb,2);
   ib_block = size(ib,2);
%
   log10KA = 12;
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
   fprintf('                    Semi-local check on block 1\n\n');
   fprintf('||Q''Q - I|| = %e', norm(Q(1:m,ilo(1):ihi(1))'*Q(1:m,ilo(1):ihi(1)) - eye(nb(1)), 'fro'));
   fprintf('                ||A - Q*R|| = %e\n\n', norm(As(1:m,1:ihi(1)) - Q(1:m,ilo(1):ihi(1))*R(ilo(1):ihi(1),ilo(1):ihi(1)), 'fro') / norm(As(1:m,1:ihi(1)), 'fro'));
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%    REMOVING LOOP TO WORK ON INDIVIDUAL BLOCKS
%    Second block
%
   lda = -1;
   ldq = -1;
   ldt = -1;
%
   ml = m;
%   for k = 2:nb_block,
%   for k = 2:ib_block,
%
      ml = ml - nb(2-1);
      nl = nb(2);
%
      [ A ] = lila_ormqrf_v02( m, nb(2), ihi(2-1), A, 1, 1, lda, A, 1, ilo(2), lda, T, 1, 1, ldt, ib, nb, 2 );
%
      [ A ] = lila_geqrf_v00( ml, nl, A, ilo(2), ilo(2), lda );
%
      T(ilo(2):ihi(2),ilo(2):ihi(2) ) = larft( A(ilo(2):m,ilo(2):ihi(2) ) );
      T(1:ihi(2-1),ilo(2):ihi(2)) = -( T(1:ihi(2-1),1:ihi(2-1))*A(ilo(2):m,1:ihi(2-1))' )*( (tril(A(ilo(2):m,ilo(2):ihi(2)),-1)+eye(m-ilo(2)+1,nb(2)))*T(ilo(2):ihi(2),ilo(2):ihi(2) ) );
%
      Q(ilo(2):m,ilo(2):ihi(2)) = A(ilo(2):m,ilo(2):ihi(2));
      [ Q ] = lila_orgqr_v02( ml, nl, Q, ilo(2), ilo(2), ldq, T, ilo(2), ilo(2), ldt );
%
      [ Q ] = lila_ormqrbz_v02( m, nb(2), ihi(2-1), A, 1, 1, lda, Q, 1, ilo(2), ldq, T, 1, 1, ldt );
%
      R = triu(A(1:ihi(2),1:ihi(2)));
      fprintf('||Q''Q - I|| = %e', norm(Q(1:m,1:ihi(2))'*Q(1:m,1:ihi(2)) - eye(ihi(2)), 'fro'));
      fprintf('                ||A - Q*R|| = %e\n\n', norm(As(1:m,1:ihi(2)) - Q(1:m,1:ihi(2))*R(1:ihi(2),1:ihi(2)), 'fro') / norm(As(1:m,1:ihi(2)), 'fro'));
return
%   end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Third block
%

      ml = ml - nb(3-1);
      nl = nb(3);
%
      [ A ] = lila_ormqrf_v02( m, nb(3), ihi(3-1), A, 1, 1, lda, A, 1, ilo(3), lda, T, 1, 1, ldt, ib, nb, 3 );
%
      [ A ] = lila_geqrf_v00( ml, nl, A, ilo(3), ilo(3), lda );
%
      T(ilo(3):ihi(3),ilo(3):ihi(3) ) = larft( A(ilo(3):m,ilo(3):ihi(3) ) );
      T(1:ihi(3-1),ilo(3):ihi(3)) = -( T(1:ihi(3-1),1:ihi(3-1))*A(ilo(3):m,1:ihi(3-1))' )*( (tril(A(ilo(3):m,ilo(3):ihi(3)),-1)+eye(m-ilo(3)+1,nb(3)))*T(ilo(3):ihi(3),ilo(3):ihi(3) ) );
%
      Q(ilo(3):m,ilo(3):ihi(3)) = A(ilo(3):m,ilo(3):ihi(3));
      [ Q ] = lila_orgqr_v02( ml, nl, Q, ilo(3), ilo(3), ldq, T, ilo(3), ilo(3), ldt );
%
      [ Q ] = lila_ormqrbz_v02( m, nb(3), ihi(3-1), A, 1, 1, lda, Q, 1, ilo(3), ldq, T, 1, 1, ldt );
%
      R = triu(A(1:ihi(3),1:ihi(3)));
      fprintf('||Q''Q - I|| = %e', norm(Q(1:m,1:ihi(3))'*Q(1:m,1:ihi(3)) - eye(ihi(3)), 'fro'));
      fprintf('                ||A - Q*R|| = %e\n\n', norm(As(1:m,1:ihi(3)) - Q(1:m,1:ihi(3))*R(1:ihi(3),1:ihi(3)), 'fro') / norm(As(1:m,1:ihi(3)), 'fro'));

%   end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Fourth block
%
      ml = ml - nb(4-1);
      nl = nb(4);
%
      [ A ] = lila_ormqrf_v02( m, nb(4), ihi(4-1), A, 1, 1, lda, A, 1, ilo(4), lda, T, 1, 1, ldt, ib, nb, 4 );
%
      [ A ] = lila_geqrf_v00( ml, nl, A, ilo(4), ilo(4), lda );
%
      T(ilo(4):ihi(4),ilo(4):ihi(4) ) = larft( A(ilo(4):m,ilo(4):ihi(4) ) );
      T(1:ihi(4-1),ilo(4):ihi(4)) = -( T(1:ihi(4-1),1:ihi(4-1))*A(ilo(4):m,1:ihi(4-1))' )*( (tril(A(ilo(4):m,ilo(4):ihi(4)),-1)+eye(m-ilo(4)+1,nb(4)))*T(ilo(4):ihi(4),ilo(4):ihi(4) ) );
%
      Q(ilo(4):m,ilo(4):ihi(4)) = A(ilo(4):m,ilo(4):ihi(4));
      [ Q ] = lila_orgqr_v02( ml, nl, Q, ilo(4), ilo(4), ldq, T, ilo(4), ilo(4), ldt );
%
      [ Q ] = lila_ormqrbz_v02( m, nb(4), ihi(4-1), A, 1, 1, lda, Q, 1, ilo(4), ldq, T, 1, 1, ldt );
%
      R = triu(A(1:ihi(4),1:ihi(4)));
      fprintf('||Q''Q - I|| = %e', norm(Q(1:m,1:ihi(4))'*Q(1:m,1:ihi(4)) - eye(ihi(4)), 'fro'));
      fprintf('                ||A - Q*R|| = %e\n\n', norm(As(1:m,1:ihi(4)) - Q(1:m,1:ihi(4))*R(1:ihi(4),1:ihi(4)), 'fro') / norm(As(1:m,1:ihi(4)), 'fro'));

%   end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Final block
%

      ml = ml - nb(5-1);
      nl = nb(5);
%
      [ A ] = lila_ormqrf_v02( m, nb(5), ihi(5-1), A, 1, 1, lda, A, 1, ilo(5), lda, T, 1, 1, ldt, ib, nb, 5 );
%
      [ A ] = lila_geqrf_v00( ml, nl, A, ilo(5), ilo(5), lda );
%
      T(ilo(5):ihi(5),ilo(5):ihi(5) ) = larft( A(ilo(5):m,ilo(5):ihi(5) ) );
      T(1:ihi(5-1),ilo(5):ihi(5)) = -( T(1:ihi(5-1),1:ihi(5-1))*A(ilo(5):m,1:ihi(5-1))' )*( (tril(A(ilo(5):m,ilo(5):ihi(5)),-1)+eye(m-ilo(5)+1,nb(5)))*T(ilo(5):ihi(5),ilo(5):ihi(5) ) );
%
      Q(ilo(5):m,ilo(5):ihi(5)) = A(ilo(5):m,ilo(5):ihi(5));
      [ Q ] = lila_orgqr_v02( ml, nl, Q, ilo(5), ilo(5), ldq, T, ilo(5), ilo(5), ldt );
%
      [ Q ] = lila_ormqrbz_v02( m, nb(5), ihi(5-1), A, 1, 1, lda, Q, 1, ilo(5), ldq, T, 1, 1, ldt );
%
      R = triu(A(1:ihi(5),1:ihi(5)));
      fprintf('||Q''Q - I|| = %e', norm(Q(1:m,1:ihi(5))'*Q(1:m,1:ihi(5)) - eye(ihi(5)), 'fro'));
      fprintf('                ||A - Q*R|| = %e\n\n', norm(As(1:m,1:ihi(5)) - Q(1:m,1:ihi(5))*R(1:ihi(5),1:ihi(5)), 'fro') / norm(As(1:m,1:ihi(5)), 'fro'));

%   end
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
