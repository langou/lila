%
   clear
%
   m = 30;
   mt = 4;
   nb = [ 25, 4 ];
%
   n = sum(nb);
   tb = mt*ones(1,n);
   if ( sum(tb) < n ) fprintf('sum(tb) < n\n'); return; end
%
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
%  Q = zeros( m, n );
   Q = randn( m, n );
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   fprintf('\n')

   TT = zeros(mt,n);

   if ( 1 + mt > nb(1) ) vt = nb(1) - 1; else vt = mt; end;

   ilot = 1;
   ihit = vt;
   jlot = 1;
   jhit = vt;

   while ( nb(1) ~= jhit )
%   fprintf('jlot=%d, jhit=%d, vt=%d\n',jlot,jhit,vt);
   TT(ilot:ihit,jlot:jhit) = T(jlot:jhit,jlot:jhit);
   if ( jhit + mt > nb(1) ) vt = nb(1) - jhit; else vt = mt; end;
   ilot = 1;
   ihit = vt;
   jlot = jhit+1;
   jhit = jhit+vt;
   end
%   fprintf('jlot=%d, jhit=%d, vt=%d\n',jlot,jhit,vt);
   TT(ilot:ihit,jlot:jhit) = T(jlot:jhit,jlot:jhit);
   if ( jhit + mt > nb(1) ) vt = nb(1) - jhit; else vt = mt; end;

%
%  Second block
%
   lda = -1;
   ldq = -1;
   ldt = -1;
%
   ml = m;
%
   jblo = ilo(2);
   jbhi = ihi(2);
   iblo = 1;
   ibhi = m;
%
   ialo = 1;
   jalo = 1;
   jahi = 1;
   iahi = m;
%
   itlo = 1;
   jtlo = 1;
   ithi = 1;
   jthi = 1;
%
   mll = m;
%
   i = 3;
   k = 10;
% 
   for ii = 1:i-1,
% 
      V = tril(A(ialo:iahi,jalo:jahi), -1) + eye(mll,1);
%
      H = (eye(mll,mll) - V * ( larft(V) * V' ) );
%
      A(iblo:ibhi,jblo:jbhi) = H*A(iblo:ibhi,jblo:jbhi);
%       
      mll = mll-1;
% 
      iblo = iblo+1;
% 
      ialo = ialo+1;
      jalo = jalo+1;
      jahi = jahi+1;
      iahi = m;
% 
      itlo = itlo+1;
      jtlo = jtlo+1;
      ithi = ithi+1;
      jthi = jthi+1;
% 
   end
%
   [ A ] = lila_ormqrf_v05( m, nb(2), k, i, A, 1, 1, lda, A, 1, ilo(2), lda, TT, mt );
%
   mll = mll-k;
%
   iblo = iblo+k;
%
   ialo = ialo+k;
   jalo = jalo+k;
   jahi = jahi+k;
   iahi = m;
%
   itlo = itlo+k;
   jtlo = jtlo+k;
   ithi = ithi+k;
   jthi = jthi+k;
%
   for ii = i+k:ihi(1),
%
      V = tril(A(ialo:iahi,jalo:jahi), -1) + eye(mll,1);
%
      H = (eye(mll,mll) - V * ( larft(V) * V' ) );
% 
      A(iblo:ibhi,jblo:jbhi) = H*A(iblo:ibhi,jblo:jbhi);
%        
      mll = mll-1;
%
      iblo = iblo+1;
%
      ialo = ialo+1;
      jalo = jalo+1;
      jahi = jahi+1;
      iahi = m;
%
      itlo = itlo+1;
      jtlo = jtlo+1;
      ithi = ithi+1;
      jthi = jthi+1;
%
   end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
   ml = ml - nb(1);
   nl = nb(2);
%
   [ A ] = lila_geqrf_v00( ml, nl, A, ilo(2), ilo(2), lda );
%
   T(ilo(2):ihi(2),ilo(2):ihi(2) ) = larft( A(ilo(2):m,ilo(2):ihi(2) ) );
   T(1:ihi(2-1),ilo(2):ihi(2)) = -( T(1:ihi(2-1),1:ihi(2-1))*A(ilo(2):m,1:ihi(2-1))' )*( (tril(A(ilo(2):m,ilo(2):ihi(2)),-1)+eye(m-ilo(2)+1,nb(2)))*T(ilo(2):ihi(2),ilo(2):ihi(2) ) );
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
   Q(ilo(2):m,ilo(2):ihi(2)) = A(ilo(2):m,ilo(2):ihi(2));
%  [ Q ] = lila_orgqr_v04( m, nb(2), ilo(2), Q, 1, 1, ldq, T, 1, 1, ldt, tb );
   V=tril(A(ilo(2):m,ilo(2):ihi(2)),-1)+eye(ml,nb(2));
   Hl=(eye(ml)-V*T(ilo(2):ihi(2),ilo(2):ihi(2))*V');
   Q(ilo(2):m,ilo(2):ihi(2)) = Hl(1:ml,1:nb(2));
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
   i = 1;
   k = 25;
%
      Q(i+k:ihi(1),ilo(2):ihi(2)) = zeros( size( Q(i+k:ihi(1),ilo(2):ihi(2)) ));
%
      ialo = ihi(1);   
      iahi = m;   
%
      jalo = ihi(1);
      jahi = ihi(1);
%
      iblo = ihi(1);
      ibhi = m;
%
      jblo = ilo(2);
      jbhi = ihi(2);
%
      itlo = ihi(1);
      ithi = ihi(1);
%
      jtlo = ihi(1);
      jthi = ihi(1);
% 
    for ii = ihi(1):-1:i+k,
% 
       ml = iahi-ialo+1;
% 
       V = tril(A(ialo:iahi,jalo:jahi), -1) + eye(size(A(ialo:iahi,jalo:jahi)));
% 
       H = (eye(ml,ml) - V * ( T(itlo:ithi,jtlo:jthi) * V' ) );
% 
       Q(iblo:ibhi,jblo:jbhi) = H*Q(iblo:ibhi,jblo:jbhi);
% 
       iblo = iblo-1;
% 
       ialo = ialo-1;
       jalo = jalo-1;
       jahi = jahi-1;
% 
       itlo = itlo-1;
       ithi = ithi-1;
% 
       jtlo = jtlo-1;
       jthi = jthi-1;
% 
    end
%
%   Q(i:i+k-1,ilo(2):ihi(2)) = zeros( size( Q(i:i+k-1,ilo(2):ihi(2)) ));
%
%if (1==0)
   [ Q ] = lila_ormqrbz_v04( m, nb(2), k, i, A, 1, 1, lda, Q, 1, ilo(2), ldq, T, 1, 1, ldt, tb );
%
   iblo = iblo-k;
%
   ialo = ialo-k;
   jalo = jalo-k;
   jahi = jahi-k;
%
   itlo = itlo-k;
   ithi = ithi-k;
%
   jtlo = jtlo-k;
   jthi = jthi-k;
%
   Q(1:i-1,ilo(2):ihi(2)) = zeros( size( Q(1:i-1,ilo(2):ihi(2)) ));
%
  for ii = i-1:-1:1,
%
      ml = iahi-ialo+1;
%
      V = tril(A(ialo:iahi,jalo:jahi), -1) + eye(size(A(ialo:iahi,jalo:jahi)));
%
      H = (eye(ml,ml) - V * ( T(itlo:ithi,jtlo:jthi) * V' ) );
%
      Q(iblo:ibhi,jblo:jbhi) = H*Q(iblo:ibhi,jblo:jbhi);
%
      iblo = iblo-1;
%
      ialo = ialo-1;
      jalo = jalo-1;
      jahi = jahi-1;
%
      itlo = itlo-1;
      ithi = ithi-1;
%
      jtlo = jtlo-1;
      jthi = jthi-1;
%
   end
%end
%
%
   R = triu(A(1:ihi(2),1:ihi(2)));
   fprintf('||Q''Q - I|| = %e', norm(Q(1:m,1:ihi(2))'*Q(1:m,1:ihi(2)) - eye(ihi(2)), 'fro'));
   fprintf('                ||A - Q*R|| / || A || = %e\n', norm(As(1:m,1:ihi(2)) - Q(1:m,1:ihi(2))*R(1:ihi(2),1:ihi(2)), 'fro') / norm(As(1:m,1:ihi(2)), 'fro'));
%
%
   fprintf('\n')
%
