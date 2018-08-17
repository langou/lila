%
   fprintf('\n');
   clear
%
   nb = [ 7, 4, 6, 10 ];
   m = 40;
%
   n = sum(nb);
   nb_block = size(nb,2);
%
   log10KA = 10;
%
   U = randn(m,n); [U,~]=qr(U,0);
   V = randn(n,n); [V,~]=qr(V,0);
   S = diag( 10.^( linspace( 0, log10KA, n ) ) );
   A = U * S * V';
   clear U S V;
%
   ilo = zeros(nb_block,1);
   ilo(1) = [1];
   for i = 2:nb_block,
      ilo(i) = ilo(i-1) + nb(i-1);
   end
%
   ihi = [nb(1)];
   for i = 2:nb_block,
      ihi(i) = ihi(i-1) + nb(i);
   end
%
   nrmA = norm(A,'fro');
   R = zeros(n);
   Q = zeros(m,n);
%
%  block gram-schmidt, PN
%
   for k = 1:4,

   Q(1:m,ilo(k):ihi(k)) = A(1:m,ilo(k):ihi(k)) ;

   if (k > 1 ), R(ilo(1):ihi(k-1),ilo(k):ihi(k)) = Q(1:m,ilo(1):ihi(k-1))'*Q(1:m,ilo(k):ihi(k)); end;

   if (k > 1 ), Q(1:m,ilo(k):ihi(k)) = Q(1:m,ilo(k):ihi(k)) - Q(1:m,ilo(1):ihi(k-1))*R(ilo(1):ihi(k-1),ilo(k):ihi(k)); end;

   [Q(1:m,ilo(k):ihi(k)) R(ilo(k):ihi(k),ilo(k):ihi(k))] = qr(Q(1:m,ilo(k):ihi(k)),0);

   end

   fprintf('|| A - Q*R || / || A ||   = %6.1e\n',norm(A(1:m,1:n)-Q(1:m,1:n)*R(1:n,1:n),'fro')/norm(A(1:m,1:n),'fro'));
   fprintf('|| I - Q''*Q ||            = %6.1e\n',norm(eye(n) - Q(1:m,1:n)'*Q(1:m,1:n),'fro'));

%
%  block gram-schmidt PNPN
%
   for k = 1:4,

   Q(1:m,ilo(k):ihi(k)) = A(1:m,ilo(k):ihi(k)) ;

   if (k > 1 ), R1(ilo(1):ihi(k-1),ilo(k):ihi(k)) = Q(1:m,ilo(1):ihi(k-1))'*Q(1:m,ilo(k):ihi(k)); end;

   if (k > 1 ), Q(1:m,ilo(k):ihi(k)) = Q(1:m,ilo(k):ihi(k)) - Q(1:m,ilo(1):ihi(k-1))*R1(ilo(1):ihi(k-1),ilo(k):ihi(k)); end;

   [Q(1:m,ilo(k):ihi(k)) R1(ilo(k):ihi(k),ilo(k):ihi(k))] = qr(Q(1:m,ilo(k):ihi(k)),0);

   if (k > 1 ), R2(ilo(1):ihi(k-1),ilo(k):ihi(k)) = Q(1:m,ilo(1):ihi(k-1))'*Q(1:m,ilo(k):ihi(k)); end;

   if (k > 1 ), Q(1:m,ilo(k):ihi(k)) = Q(1:m,ilo(k):ihi(k)) - Q(1:m,ilo(1):ihi(k-1))*R2(ilo(1):ihi(k-1),ilo(k):ihi(k)); end;

   [Q(1:m,ilo(k):ihi(k)) R2(ilo(k):ihi(k),ilo(k):ihi(k))] = qr(Q(1:m,ilo(k):ihi(k)),0);

   if (k > 1 ), R(ilo(1):ihi(k-1),ilo(k):ihi(k)) = R1(ilo(1):ihi(k-1),ilo(k):ihi(k)) + R2(ilo(1):ihi(k-1),ilo(k):ihi(k)) * R1(ilo(k):ihi(k),ilo(k):ihi(k)); end;

   if (k > 1 ), R(ilo(k):ihi(k),ilo(k):ihi(k)) = R2(ilo(k):ihi(k),ilo(k):ihi(k)) * R1(ilo(k):ihi(k),ilo(k):ihi(k));  end;

   end

   fprintf('|| A - Q*R || / || A ||   = %6.1e\n',norm(A(1:m,1:n)-Q(1:m,1:n)*R(1:n,1:n),'fro')/norm(A(1:m,1:n),'fro'));
   fprintf('|| I - Q''*Q ||            = %6.1e\n',norm(eye(n) - Q(1:m,1:n)'*Q(1:m,1:n),'fro'));

%
%  block gram-schmidt PPN
%
   for k = 1:4,

   Q(1:m,ilo(k):ihi(k)) = A(1:m,ilo(k):ihi(k)) ;

   if (k > 1 ), R1(ilo(1):ihi(k-1),ilo(k):ihi(k)) = Q(1:m,ilo(1):ihi(k-1))'*Q(1:m,ilo(k):ihi(k)); end;

   if (k > 1 ), Q(1:m,ilo(k):ihi(k)) = Q(1:m,ilo(k):ihi(k)) - Q(1:m,ilo(1):ihi(k-1))*R1(ilo(1):ihi(k-1),ilo(k):ihi(k)); end;

   if (k > 1 ), R2(ilo(1):ihi(k-1),ilo(k):ihi(k)) = Q(1:m,ilo(1):ihi(k-1))'*Q(1:m,ilo(k):ihi(k)); end;

   if (k > 1 ), Q(1:m,ilo(k):ihi(k)) = Q(1:m,ilo(k):ihi(k)) - Q(1:m,ilo(1):ihi(k-1))*R2(ilo(1):ihi(k-1),ilo(k):ihi(k)); end;

   [Q(1:m,ilo(k):ihi(k)) R(ilo(k):ihi(k),ilo(k):ihi(k))] = qr(Q(1:m,ilo(k):ihi(k)),0);

   if (k > 1 ), R(ilo(1):ihi(k-1),ilo(k):ihi(k)) = R1(ilo(1):ihi(k-1),ilo(k):ihi(k)) + R2(ilo(1):ihi(k-1),ilo(k):ihi(k)); end;

   end

   fprintf('|| A - Q*R || / || A ||   = %6.1e\n',norm(A(1:m,1:n)-Q(1:m,1:n)*R(1:n,1:n),'fro')/norm(A(1:m,1:n),'fro'));
   fprintf('|| I - Q''*Q ||            = %6.1e\n',norm(eye(n) - Q(1:m,1:n)'*Q(1:m,1:n),'fro'));












return











%   R1 = chol( A(1:m,1:nb(1))'*A(1:m,1:nb(1)), 'upper'); Q1 = A(1:m,1:nb(1)) / R1;
%  R1 = chol( A(1:m,ilo(1):ihi(1))'*A(1:m,ilo(1):ihi(1)), 'upper'); Q1 = A(1:m,ilo(1):ihi(1)) / R1;
%
   fprintf('block one is of size %d', m)
   fprintf('x%d', nb(1));
   fprintf('\n'); 
   fprintf('|| A1 - Q1*R1 || / || A ||   = %6.1e\n',norm(A(1:m,ilo(1):ihi(1))-Q1*R1,'fro')/norm(A(1:m,ilo(1):ihi(1)),'fro'));
   fprintf('|| I1 - Q1''*Q1 ||            = %6.1e\n',norm(eye(ihi(1)) - Q1'*Q1,'fro'));
%
   for i=ilo(1):ihi(1),
      if( R1(i,i) < 0 ), Q1(:,i) = - Q1(:,i); R1(i,:) = - R1(i,:); end
      if( rand(1) >  0.5 ), Q1(:,i) = - Q1(:,i); R1(i,:) = - R1(i,:); end
   end
%
   Qs_1 = Q1; Rs_1 = R1;
%
   D1 = -eye(ihi(1),ihi(1));
   for k = ilo(1):ihi(1),
      if (abs(1 - Q1(k,k)) < abs( 1 + Q1(k,k) )) 
      else 
         Q1(1:m,k) = -Q1(1:m,k);
         D1(k,k) = 1.0e+00;
      end
      Q1(k,k) = 1 + Q1(k,k);
      Q1(k+1:m,k) = Q1(k+1:m,k) / Q1(k,k);
      Q1(k+1:m,k+1:ihi(1)) = Q1(k+1:m,k+1:ihi(1)) -  Q1(k+1:m,k) * Q1(k,k+1:ihi(1));
   end
%
%  Checks and construction of H
%
   U1 = triu(Q1); U1=U1(1:ihi(1),1:ihi(1));
   L1 = tril(Q1,-1)+eye(m,ihi(1));
   fprintf('|| Q1 - L1*U1 ||             = %d \n', norm( ( eye(m,ihi(1)) - Qs_1 * D1 ) - L1*U1,'fro') );
%
   V1 = L1;
   T1 = U1 / (V1(1:ihi(1),1:ihi(1))'); %Not used currently
   TT1 = larft(V1);
%
   H1 = ( eye(m) - V1 * TT1 * V1' );
%
   fprintf('\n');
   fprintf('||H1''*H1 - I|| = %e (1,1) \n', norm(H1(1:m,1:ihi(1))'*H1(1:m,1:ihi(1)) - eye(ihi(1)),'fro'));
   fprintf('||H1''*H1 - I|| = %e (2,2) \n', norm(H1(1:m,m-ihi(1):m)'*H1(1:m,m-ihi(1):m) - eye(m-ihi(1)),'fro'));
   fprintf('||H1''*H1 - I|| = %e (1,2) \n', norm(H1(1:m,1:ihi(1))'*H1(1:m,m-ihi(1):m) ,'fro'));
   fprintf('||H1''*H1 - I|| = %e (all) \n\n', norm(H1'*H1 - eye(m),'fro'));

   fprintf(' ----- ||H1''*A - R|| / ||R|| = %e\n', ( norm( triu(H1(1:m,1:ihi(1))'*A(1:m,ihi(1))) - D1*Rs_1,'fro') / nrmA ) );
   fprintf(' ----- ||tril(H1''*A, -1)||   = %e\n', ( norm(tril(H1'*A(1:m,ihi(1)), -1),'fro' ) )  / nrmA );
   fprintf(' ----- ||H1 - Q||            = %e', ( norm(H1(1:m,1:ihi(1)) - Qs_1*D1,'fro') ) );
   fprintf('\n\n');
   fprintf('||T - larft(V)|| = %e ', norm(T1 - TT1,'fro'));
   fprintf('\n\n');
%







return
   n1 = 4;
   B = randn(m,n1);
%
   B1 = B - Qs*(Qs'*B);
   fprintf('norm( Qs'' * B1 ) / norm( B1 ) = %8.1e\n', norm( Qs' * B1 ) / norm( B1 ));
   [GS_B1] = geqr2(B1);
   TT_GS_B1 = larft(GS_B1);
   V_GS_B1 = tril(GS_B1,-1) + eye( m, n1);

   ( eye(m,m) - V_GS_B1 * TT_GS_B1' * V_GS_B1' ) * B1;

   Q_GS_2 = ( eye(m,m) - V_GS_B1 * TT_GS_B1 * V_GS_B1' ) * eye(m,n1);
  
%
%
%
%  B11 = ( B - V * ( TT * ( V' * B ) ) );
%  B22 = ( B - V * ( T * ( V' * B ) ) );
%
%  norm( H(1:m,1:n1)' * B11 ) / norm( B11 ), norm( H(1:m,1:n1)' * B22 ) / norm( B22 )
%  norm( Qs' * B11 ) / norm( B11 ), norm( Qs' * B22 ) / norm( B22 )
%
   HB = H1'*B;
   HB_1 = HB(1:n,1:n1);
   AA_2 = HB(n+1:m,1:n1);
   HB_2 = HB(n+1:m,1:n1);

   [HB_2] = geqr2(HB_2)
   TT2 = larft(HB_2);
   V_2 = tril(HB_2,-1) + eye( m-n, n1);

   ( eye(m-n,m-n) - V_2 * TT2' * V_2' ) * AA_2
   Q_HH_2_small = ( eye(m-n,m-n) - V_2 * TT2 * V_2' ) * eye(m-n,n1);

   Q_HH_2 = H1 * [ zeros(n,n1); Q_HH_2_small ];

   





   fprintf('\n');
   diag(D)';
