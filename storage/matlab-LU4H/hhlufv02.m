%
   fprintf('\n');
   clear
%
   m = 15;
   n = 7;
%
   log10KA = 1.5;
%
   U = randn(m,n); [U,~]=qr(U,0);
   V = randn(n,n); [V,~]=qr(V,0);
   S = diag( 10.^( linspace( 0, log10KA, n ) ) );
   A = U * S * V';
   clear U S V;
%
   nrmA = norm(A,'fro');
%
  [Q R] = qr(A,0);
%   R = chol( A'*A, 'upper'); Q = A / R;
%
   fprintf('|| A - Q*R || / || A ||   = %6.1e\n',norm(A-Q*R,'fro')/norm(A,'fro'));
   fprintf('|| I - Q''*Q ||            = %6.1e\n',norm(eye(n) - Q'*Q,'fro'));
%
   for i=1:n,
      if( R(i,i) < 0 ), Q(:,i) = - Q(:,i); R(i,:) = - R(i,:); end
      if( rand(1) >  0.5 ), Q(:,i) = - Q(:,i); R(i,:) = - R(i,:); end
   end
%
   Qs = Q; Rs = R;
%
   D = -eye(n,n);
   for k = 1:n,
      if (abs(1 - Q(k,k)) < abs( 1 + Q(k,k) )) 
      else 
         Q(1:m,k) = -Q(1:m,k);
         D(k,k) = 1.0e+00;
      end
      Q(k,k) = 1 + Q(k,k);
      Q(k+1:m,k) = Q(k+1:m,k) / Q(k,k);
      Q(k+1:m,k+1:n) = Q(k+1:m,k+1:n) -  Q(k+1:m,k) * Q(k,k+1:n);
   end
%
%  Checks and construction of H
%
   U = triu(Q); U=U(1:n,1:n);
   L = tril(Q,-1)+eye(m,n);
   fprintf('|| Q - L*U ||             = %d \n', norm( ( eye(m,n) - Qs * D ) - L*U,'fro') );
%
   V = L;
   T = U / (V(1:n,1:n)');
   TT = larft(V); T = TT;
%
   H = ( eye(m) - V * T * V' );
   H1 = ( eye(m) - V * TT * V' );
%
   fprintf('\n');
   fprintf('||H''*H - I||   = %e (1,1) \n', norm(H(1:m,1:n)'*H(1:m,1:n) - eye(n),'fro'));
   fprintf('||H''*H - I||   = %e (2,2) \n', norm(H(1:m,m-n:m)'*H(1:m,m-n:m) - eye(m-n),'fro'));
   fprintf('||H''*H - I||   = %e (1,2) \n', norm(H(1:m,1:n)'*H(1:m,m-n:m) ,'fro'));
   fprintf('||H''*H - I||   = %e (all) \n', norm(H'*H - eye(m),'fro'));
   fprintf('||H1''*H1 - I|| = %e (1,1) \n', norm(H1(1:m,1:n)'*H1(1:m,1:n) - eye(n),'fro'));
   fprintf('||H1''*H1 - I|| = %e (2,2) \n', norm(H1(1:m,m-n:m)'*H1(1:m,m-n:m) - eye(m-n),'fro'));
   fprintf('||H1''*H1 - I|| = %e (1,2) \n', norm(H1(1:m,1:n)'*H1(1:m,m-n:m) ,'fro'));
   fprintf('||H1''*H1 - I|| = %e (all) \n\n', norm(H1'*H1 - eye(m),'fro'));

   fprintf(' ----- ||H''*A - R|| / ||R||  = %e\n', ( norm( triu(H(1:m,1:n)'*A) - D*Rs,'fro') / nrmA ) );
   fprintf(' ----- ||tril(H''*A, -1)||    = %e\n', ( norm(tril(H'*A, -1),'fro' ) )  / nrmA );
   fprintf(' ----- ||H - Q||             = %e \n', ( norm(H(1:m,1:n) - Qs*D,'fro') ) );
   fprintf(' ----- ||H1''*A - R|| / ||R|| = %e\n', ( norm( triu(H1(1:m,1:n)'*A) - D*Rs,'fro') / nrmA ) );
   fprintf(' ----- ||tril(H1''*A, -1)||   = %e\n', ( norm(tril(H1'*A, -1),'fro' ) )  / nrmA );
   fprintf(' ----- ||H1 - Q||            = %e', ( norm(H1(1:m,1:n) - Qs*D,'fro') ) );
   fprintf('\n\n');
   fprintf('||T - larft(V)|| = %e ', norm(T - TT,'fro'));
   fprintf('\n\n');
%
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
