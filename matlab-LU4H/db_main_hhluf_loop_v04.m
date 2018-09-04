%
   clear
%
   m = 60;
  tb = [ 4, 5, 3, 5, 7, 9, 5, 4, 5, 7 ];
  nb = [ 20, 13 ];
%   tb = [ 10, 10:-2:1 ];
%   nb = [ 10, 35 ];
   n = sum(nb);
   if ( sum(tb) < n ) fprintf('sum(tb) < n\n'); return; end
% tb = [ 4 ];
% for i = 1:sum(tb),
% for k = 1:sum(tb)-i+1,
%   tb = 50;
%
   for i_1 = 3,
   for k_1 = 6,
%
   nb_block = size(nb,2);
   tb_block = size(tb,2);
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
   for ii = 2:nb_block,
      ilo(ii) = ilo(ii-1) + nb(ii-1);
   end
%
   ihi = zeros(nb_block,1);
   ihi = nb(1);
   for ii = 2:nb_block,
      ihi(ii) = ihi(ii-1) + nb(ii);
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
%  fprintf('Local check on block 1');
%  fprintf('                    Semi-local check on block 1\n\n');
%  fprintf('||Q''Q - I|| = %e', norm(Q(1:m,ilo(1):ihi(1))'*Q(1:m,ilo(1):ihi(1)) - eye(nb(1)), 'fro'));
%  fprintf('                ||A - Q*R|| = %e\n\n', norm(As(1:m,1:ihi(1)) - Q(1:m,ilo(1):ihi(1))*R(ilo(1):ihi(1),ilo(1):ihi(1)), 'fro') / norm(As(1:m,1:ihi(1)), 'fro'));
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   ml = m;
%

   jblo = ilo(2);
   jbhi = ihi(2);
   iblo = 1;
   ibhi = m;

   ialo = 1;
   jalo = 1;
   jahi = 1;
   iahi = m;


   itlo = 1;
   jtlo = 1;
   ithi = 1;
   jthi = 1;

   mll = m;

   i = i_1;
   k = k_1;
 
   for ii = 1:i-1,
 
      V = tril(A(ialo:iahi,jalo:jahi), -1) + eye(mll,1);

      H = (eye(mll,mll) - V * ( larft(V) * V' ) );

      A(iblo:ibhi,jblo:jbhi) = H*A(iblo:ibhi,jblo:jbhi);
       
      mll = mll-1;
 
      iblo = iblo+1;
 
      ialo = ialo+1;
      jalo = jalo+1;
      jahi = jahi+1;
      iahi = m;
 
      itlo = itlo+1;
      jtlo = jtlo+1;
      ithi = ithi+1;
      jthi = jthi+1;
 
   end

   [ A ] = lila_ormqrf_v03( m, nb(2), k, i, A, 1, 1, lda, A, 1, ilo(2), lda, T, 1, 1, ldt, tb );

   mll = mll-k;

   iblo = iblo+k;

   ialo = ialo+k;
   jalo = jalo+k;
   jahi = jahi+k;
   iahi = m;

   itlo = itlo+k;
   jtlo = jtlo+k;
   ithi = ithi+k;
   jthi = jthi+k;


   for ii = i+k:ihi(1),

      V = tril(A(ialo:iahi,jalo:jahi), -1) + eye(mll,1);

      H = (eye(mll,mll) - V * ( larft(V) * V' ) );
 
      A(iblo:ibhi,jblo:jbhi) = H*A(iblo:ibhi,jblo:jbhi);
        
      mll = mll-1;

      iblo = iblo+1;

      ialo = ialo+1;
      jalo = jalo+1;
      jahi = jahi+1;
      iahi = m;

      itlo = itlo+1;
      jtlo = jtlo+1;
      ithi = ithi+1;
      jthi = jthi+1;

   end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
   ml = ml - nb(1);
   nl = nb(2);
%
%
   [ A ] = lila_geqrf_v00( ml, nl, A, ilo(2), ilo(2), lda );
%
   T(ilo(2):ihi(2),ilo(2):ihi(2) ) = larft( A(ilo(2):m,ilo(2):ihi(2) ) );
   T(1:ihi(2-1),ilo(2):ihi(2)) = -( T(1:ihi(2-1),1:ihi(2-1))*A(ilo(2):m,1:ihi(2-1))' )*( (tril(A(ilo(2):m,ilo(2):ihi(2)),-1)+eye(m-ilo(2)+1,nb(2)))*T(ilo(2):ihi(2),ilo(2):ihi(2) ) );
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%
%               Working to break this bad boy up now
%

   Q(ilo(2):m,ilo(2):ihi(2)) = A(ilo(2):m,ilo(2):ihi(2));
   [ Q ] = lila_orgqr_v03( m, nb(2), ilo(2), Q, 1, 1, ldq, T, 1, 1, ldt, tb );


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%


%%%
%  [ Q ] = lila_ormqrbz_v02( m, nb(2), ihi(2-1), A, 1, 1, lda, Q, 1, ilo(2), ldq, T, 1, 1, ldt );
%
   i = 4;
   k = 10;
%
   Q(1:ihi(1),ilo(2):ihi(2)) = zeros( size( Q(1:ihi(1),ilo(2):ihi(2)) ));

      ialo = ihi(1);   % k = ihi(1);
%
      iahi = m;    % Stays the same
%
      jalo = ihi(1);
      jahi = ihi(1);
%
      iblo = ihi(1);
      ibhi = m;
      jblo = ilo(2);
%
      jbhi = ihi(2);
%
      itlo = ihi(1);
      ithi = ihi(1);
      jtlo = ihi(1);
      jthi = ihi(1);
%
   for ii = ihi(1):-1:i+k,

      ml = iahi-ialo+1;

      V = tril(A(ialo:iahi,jalo:jahi), -1) + eye(size(A(ialo:iahi,jalo:jahi)));

      H = (eye(ml,ml) - V * ( T(itlo:ithi,jtlo:jthi) * V' ) );

      Q(iblo:ibhi,jblo:jbhi) = H*Q(iblo:ibhi,jblo:jbhi);

      iblo = iblo-1;

      ialo = ialo-1;
      jalo = jalo-1;
      jahi = jahi-1;

      itlo = itlo-1;
      ithi = ithi-1;

      jtlo = jtlo-1;
      jthi = jthi-1;

   end

   [ Q ] = lila_ormqrbz_v03( m, nb(2), k, i, A, 1, 1, lda, Q, 1, ilo(2), ldq, T, 1, 1, ldt, tb );

   iblo = iblo-k;

   ialo = ialo-k;
   jalo = jalo-k;
   jahi = jahi-k;

   itlo = itlo-k;
   ithi = ithi-k;

   jtlo = jtlo-k;
   jthi = jthi-k;

  for ii = i-1:-1:1,

      ml = iahi-ialo+1;


      V = tril(A(ialo:iahi,jalo:jahi), -1) + eye(size(A(ialo:iahi,jalo:jahi)));

      H = (eye(ml,ml) - V * ( T(itlo:ithi,jtlo:jthi) * V' ) );


      Q(iblo:ibhi,jblo:jbhi) = H*Q(iblo:ibhi,jblo:jbhi);

      iblo = iblo-1;

      ialo = ialo-1;
      jalo = jalo-1;
      jahi = jahi-1;

      itlo = itlo-1;
      ithi = ithi-1;

      jtlo = jtlo-1;
      jthi = jthi-1;

   end






   R = triu(A(1:ihi(2),1:ihi(2)));
   fprintf('||Q''Q - I|| = %e', norm(Q(1:m,1:ihi(2))'*Q(1:m,1:ihi(2)) - eye(ihi(2)), 'fro'));
   fprintf('                ||A - Q*R|| / || A || = %e\n', norm(As(1:m,1:ihi(2)) - Q(1:m,1:ihi(2))*R(1:ihi(2),1:ihi(2)), 'fro') / norm(As(1:m,1:ihi(2)), 'fro'));



   end
   end



   fprintf('\n')



return

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
