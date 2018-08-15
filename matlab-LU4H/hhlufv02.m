%
   fprintf('\n');
   clear
%
   m = 15;
   n = 7;
%
   log10KA = 2;
%
   U = randn(m,n); [U,~]=qr(U,0);
   V = randn(n,n); [V,~]=qr(V,0);
   S = diag( 10.^( linspace( 0, log10KA, n ) ) );
   A = U * S * V';
   clear U S V;
%
   nrmA = norm(A,'fro');
%
%  [Q R] = qr(A,0);
   R = chol( A'*A, 'upper'); Q = A / R;
%
   fprintf('|| A - Q*R || / || A ||          = %6.1e\n',norm(A-Q*R,'fro')/norm(A,'fro'));
   fprintf('|| I - Q''*Q ||                   = %6.1e\n',norm(eye(n) - Q'*Q,'fro'));
%
   for i=1:n,
      if( R(i,i) < 0 ), Q(:,i) = - Q(:,i); R(i,:) = - R(i,:); end
      if( rand(1) >  0.5 ), Q(:,i) = - Q(:,i); R(i,:) = - R(i,:); end
   end
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
   fprintf('|| Q - L*U ||  = %d', norm( ( eye(m,n) - Qs * D ) - L*U,'fro') );
%
   V = L;
   T = U / (V(1:n,1:n)');
   TT = larft(V);
%
   H = ( eye(m) - V * T * V' );
%
   fprintf('\n');
   fprintf('||H''*H - I|| = %d (1,1) \n', norm(H(1:m,1:n)'*H(1:m,1:n) - eye(n),'fro'));
   fprintf('||H''*H - I|| = %d (2,2) \n', norm(H(1:m,m-n:m)'*H(1:m,m-n:m) - eye(m-n),'fro'));
   fprintf('||H''*H - I|| = %d (1,2) \n', norm(H(1:m,1:n)'*H(1:m,m-n:m) ,'fro'));
   fprintf('||H''*H - I|| = %d (all) \n', norm(H'*H - eye(m),'fro'));
   fprintf(' ----- ||H''*A - R|| / ||R|| = %d\n', ( norm( triu(H(1:m,1:n)'*A) - D*Rs,'fro') / nrmA ) );
   fprintf(' ----- ||tril(H''*A, -1)|| = %d\n', ( norm(tril(H'*A, -1),'fro' ) )  / nrmA );
   fprintf(' ----- ||H - Q|| = %d', ( norm(H(1:m,1:n) - Qs*D,'fro') ) );
   fprintf('\n');
   fprintf('||T - larft(V)|| = %d ', norm(T - TT,'fro'));
%
   n1 = 20;
   B = randn(m,n1);
   B1 = B - Qs*(Qs'*B);
   norm( Qs' * B1 ) / norm( B1 )
   B1 = ( B - V * ( T * ( V' * B ) ) );
   norm( Qs' * B1 ) / norm( B1 )
   B1 = ( B - V * ( TT * ( V' * B ) ) );
   norm( Qs' * B1 ) / norm( B1 )
   




   fprintf('\n');
   diag(D)'
