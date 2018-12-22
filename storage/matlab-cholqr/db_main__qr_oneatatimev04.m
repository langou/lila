%
   clear
%
   m = 76;
   n = 27;
   log10KA = 4;
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
   [ RR ] = qr(A(1:m,1:n));
%
%  all quantites starting with a d must be computed with A only: dtau, dalpha, dbeta, dgamma, etc.
%
   dA(1:n,1:n) = A(1:n,1:n);
   dV(1:n,1:n) = dA(1:n,1:n);
   dGG(1:n,1:n) = A(1:m,1:n)' * A(1:m,1:n);
%
   for k = 1:n,
%
   GG(k:n,k:n) = V(k:m,k:n)'*V(k:m,k:n);
%
   norma = norm(V(k:m,k),2);
   if ( V(k,k) > 0 ) V(k,k) = V(k,k) + norma; else V(k,k) = V(k,k) - norma; end
   alpha(k) = V(k,k);
   V(k+1:m,k) = V(k+1:m,k) / alpha(k);
   if ( V(k,k) > 0 ) V(k,k) = - norma; else V(k,k) = norma; end
   tau(k) = 2/( 1 + norm(V(k+1:m,k))^2 );
   V(k:m,k+1:n) = V(k:m,k+1:n) - [ [ 1; V(k+1:m,k) ] ] * tau(k) * ( [ 1; V(k+1:m,k) ]' * V(k:m,k+1:n) );
%
   fprintf('  R :: %e\n',norm( triu(RR(k,k:n)) - triu(V(k,k:n)), 'fro' ) / norm( triu(RR(k,k:n)), 'fro' ));
   fprintf('  R :: %e\n',norm( tril(RR(k:m,k),-1) - tril(V(k:m,k),-1), 'fro' ) / norm( tril(RR(k:m,k),-1), 'fro' ));
%
   if (k>1) dGG(k-1:n,k-1:n) = dGG(k-1:n,k-1:n) - dV(k-1,k-1:n)' * dV(k-1,k-1:n); end
   fprintf('  0 :: %e\n',norm ( dGG(k:n,k:n) - GG(k:n,k:n) , 'fro' ) / norm( GG(k:n,k:n), 'fro' ));
%
   dnorma = sqrt( dGG(k,k) );
   fprintf('  1 :: %e\n',abs( dnorma - norma ) / abs( norma ) );
%
   dbeta(k) = dV(k,k);
   dGG(k,k) = dGG(k,k) - dV(k,k)^2;
   if ( dV(k,k) > 0.e+00 ) dV(k,k) = dV(k,k) + norma; else dV(k,k) = dV(k,k) - norma; end
   dalpha(k) = dV(k,k);
   dGG(k,k) = dGG(k,k) + dV(k,k)^2;
%
   dV(k+1:n,k) = dV(k+1:n,k) / dalpha(k);
%
   if ( dV(k,k) > 0.e+00 ) dV(k,k) = - norma; else dV(k,k) = norma; end
%
   dtau(k) = 2.0e+00 / ( ( dGG(k,k) ) / dalpha(k)^2 );
   fprintf('  6 :: %e\n',abs( dtau(k) - tau(k) ) / abs( tau(k) ) );
%
   dV(k:n,k+1:n) = dV(k:n,k+1:n) - [ 1; dV(k+1:n,k) ] * dtau(k) * ( dV(k,k+1:n) + ( dGG(k,k+1:n) - dbeta(k)' * dV(k,k+1:n) ) / dalpha(k) );
%
   fprintf('  7 :: %e\n',norm( triu(RR(k,1:n)) - triu(dV(k,1:n)), 'fro' ) / norm( triu(RR(k,1:n)), 'fro' ));
   fprintf('  7 :: %e\n',norm( tril(RR(1:n,k),-1) - tril(dV(1:n,k),-1), 'fro' ) / norm( tril(RR(1:n,k),-1), 'fro' ));
   fprintf('  7 :: %e\n',norm( V(1:n,1:n) - dV(1:n,1:n), 'fro' ) / norm( V(1:n,1:n), 'fro' ));
   fprintf('  7 :: %e\n',norm( RR(1:n,1:n) - dV(1:n,1:n), 'fro' ) / norm( V(1:n,1:n), 'fro' ));
   end
   fprintf('\n');
   fprintf('\n');
%
%
   R = triu( V(1:n,1:n) );
   Q = As / R;
   fprintf('|| A - Q*R || / || A || = %6.1e\n',norm(As-Q*R)/norm(As));
   fprintf('|| I - Q''*Q || = %6.1e\n',norm(eye(n) - Q'*Q));
return


   
%
   %x = triu( dV(1:n-1,1:n-1) ) \ 





