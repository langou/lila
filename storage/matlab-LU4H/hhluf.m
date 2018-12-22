%
   fprintf('\n');
   clear
%
   m = 13;
   n = 7;
%
   log10KA = 10;
%
   U = randn(m,n); [U,~]=qr(U,0);
   V = randn(n,n); [V,~]=qr(V,0);
   S = diag( 10.^( linspace( 0, log10KA, n ) ) );
   A = U * S * V';
   clear U S V;
%
   A =  randn(m,n)*(50);
%
   [Q R] = qr(A,0);
   nrmA = norm(A,'fro');
   Q = -Q; R =-R;
%
   Q1 = Q(1:n,1:n);
   Q2 = Q(n+1:m,1:n);
%
   Q = eye(m,n) - Q;
   Qs = Q;
   nrmQ = norm(Q,'fro');
%
   sign__ = zeros(1,n);
   for k = 1:n,
      if (abs(Q(k,k)) < abs(Q(k,k) - 2)) 
         Q(k,k) = Q(k,k) + 2;
         sign__(k) = +2;
      end
      Q(k+1:m,k) = Q(k+1:m,k) / Q(k,k);
      Q(k+1:m,k+1:n) = Q(k+1:m,k+1:n) -  Q(k+1:m,k) * Q(k,k+1:n);
   end
%
%  Checks and construction of H
%
   U = triu(Q); U=U(1:n,1:n);
   L = tril(Q,-1)+eye(m,n);
   fprintf('|| Q - L*U ||  = %d', norm((Qs+[diag(sign__(1,1:n));zeros(m-n,n)] ) - L*U,'fro') / nrmQ );
%
   V = L;
   T = U / (V(1:n,1:n)');
%
   H = ( eye(m) - V * T * V' );
   [Q,R] = qr(A,0);
   Q = -Q; R =-R;

%
   fprintf('\n');
   fprintf('||H''*H - I|| = %d', norm(H'*H - eye(m),'fro'));
   fprintf(' ----- ||H''*A - R|| / ||R|| = %d', ( norm( triu(H(1:m,1:n)'*A) - R,'fro') / nrmA ) );
   fprintf(' ----- ||tril(H''*A, -1)|| = %d', ( norm(tril(H'*A, -1),'fro' ) )  / nrmA );
   fprintf(' ----- ||H''*Q - I|| = %d', ( norm(H'*Q - eye(m,n),'fro') ) );
   fprintf('\n');

   fprintf('\n');
