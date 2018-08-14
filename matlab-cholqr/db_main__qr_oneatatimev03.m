%
   clear
%
   m = 20;
   n = 4;
   log10KA = 2;
%
   U = randn(m,n); [U,~]=qr(U,0);
   V = randn(n,n); [V,~]=qr(V,0);
   S = diag( 10.^( linspace( 0, log10KA, n ) ) );
   A = U * S * V';
   clear U S V;
%
   Q = zeros(m,n);
   R = zeros(n,n);
   As = A;
   V = A;
%
   [ ~, RR ] = qr(As(1:m,1:n),0);
%
%  all quantites starting with a d must be computed with A only: dtau, dalpha, dbeta, dgamma, etc.
%
%  step 1
%  this is what we would do working on a V
   norma = norm(V(1:m,1),2); 
   if ( V(1,1) > 0 ) V(1,1) = V(1,1) + norma; else V(1,1) = V(1,1) - norma; end
   V(2:m,1) = V(2:m,1) / V(1,1);
   if ( V(1,1) > 0 ) V(1,1) = - norma; else V(1,1) = norma; end
   tau1 = 2/( 1 + norm(V(2:m,1))^2 );
   V(1,2:n) = A(1,2:n) - [ 1 ] * tau1 * ( [ 1; V(2:m,1) ]' * A(1:m,2:n) );
%
%  this is the check
   norm( RR(1,1:n) - triu(V(1,1:n)), 'fro' ) / norm( RR(1,1:n), 'fro' )
%
%  step 1
%  this is what we do working on A directly
   dA(1:n,1:n) = A(1:n,1:n);
   dG(1:n,1:n) = A(n+1:m,1:n)' * A(n+1:m,1:n);
%
   norma = sqrt( dA(1:n,1)'*dA(1:n,1) + dG(1,1) ); 
%  V(1:m,1) = A(1:m,1);
   if ( dA(1,1) > 0 ) dA(1,1) = dA(1,1) + norma; else dA(1,1) = dA(1,1) - norma; end
   dalpha(1) = dA(1,1);
%  V(2:m,1) = V(2:m,1) / dalpha(1);
   if ( dA(1,1) > 0 ) dA(1,1) = - norma; else dA(1,1) = norma; end
%
%  dtau(1) = 2/( 1 + norm(V(2:m,1))^2 );
   dtau(1) = 2/( 1 + ( dA(2:n,1)'*dA(2:n,1) + dG(1,1) ) / dalpha(1)^2 );
%
%  dA(1,2:n) = dA(1,2:n) - [ 1 ] * tau1 * ( [ 1; V(2:m,1) ]' * A(1:m,2:n) );
%  dA(1,2:n) = dA(1,2:n) - [ 1 ] * tau1 * ( dA(1,2:n) + [ V(2:m,1) ]' * A(2:m,2:n) );
%  dA(1,2:n) = dA(1,2:n) - [ 1 ] * tau1 * ( dA(1,2:n) + ( dG(1,2:n) + dA(2:n,1)' * dA(2:n,2:n) ) / dalpha(1) );
   dA(1,2:n) = dA(1,2:n) - dtau(1) * ( dA(1,2:n) + ( dG(1,2:n) + dA(2:n,1)' * dA(2:n,2:n) ) / dalpha(1) );
%
%  this is the check
   norm( RR(1,1:n) - triu(dA(1,1:n)), 'fro' ) / norm( RR(1,1:n), 'fro' )
