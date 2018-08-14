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
   [ RR ] = qr(As(1:m,1:n));
%
%  all quantites starting with a d must be computed with A only: dtau, dalpha, dbeta, dgamma, etc.
%
   dA(1:n,1:n) = A(1:n,1:n);
   dG(1:n,1:n) = A(n+1:m,1:n)' * A(n+1:m,1:n);
   dV(1:n,1:n) = dA(1:n,1:n);
%
%  step 1
%  this is what we would do working on a V
   norma = norm(V(1:m,1),2);                                                       % 1.1
   if ( V(1,1) > 0 ) V(1,1) = V(1,1) + norma; else V(1,1) = V(1,1) - norma; end    % 1.2
   alpha(1) = V(1,1);                                                              % 1.3
   V(2:m,1) = V(2:m,1) / alpha(1);                                                 % 1.4
   if ( V(1,1) > 0 ) V(1,1) = - norma; else V(1,1) = norma; end                    % 1.5
   tau1 = 2/( 1 + norm(V(2:m,1))^2 );                                              % 1.6
   V(1:m,2:n) = A(1:m,2:n) - [ [ 1; V(2:m,1) ] ] * tau1 * ( [ 1; V(2:m,1) ]' * A(1:m,2:n) );         % 1.7
%
%  this is the check
   norm( triu(RR(1,1:n)) - triu(V(1,1:n)), 'fro' ) / norm( triu(RR(1,1:n)), 'fro' )
   norm( tril(RR(1:m,1),-1) - tril(V(1:m,1),-1), 'fro' ) / norm( tril(RR(1:m,1),-1), 'fro' )
%
%  step 1
%  this is what we do working on A directly
%
%  norma = norm(V(1:m,1),2);                                                       % 1.1
   norma = sqrt( dA(1:n,1)'*dA(1:n,1) + dG(1,1) ); 
%
%  if ( V(1,1) > 0 ) V(1,1) = V(1,1) + norma; else V(1,1) = V(1,1) - norma; end    % 1.2
   if ( dV(1,1) > 0 ) dV(1,1) = dV(1,1) + norma; else dV(1,1) = dV(1,1) - norma; end
%
%  alpha(1) = V(1,1);                                                              % 1.3
   dalpha(1) = dV(1,1);
%
%  V(2:m,1) = V(2:m,1) / dalpha(1);                                                % 1.4
   dV(2:n,1) = dV(2:n,1) / dalpha(1);
%
%  if ( V(1,1) > 0 ) V(1,1) = - norma; else V(1,1) = norma; end                    % 1.5
   if ( dV(1,1) > 0 ) dV(1,1) = - norma; else dV(1,1) = norma; end
%
%  dtau(1) = 2/( 1 + norm(V(2:m,1))^2 );                                           % 1.6
   dtau(1) = 2/( 1 + dV(2:n,1)'*dV(2:n,1) + ( dG(1,1) ) / dalpha(1)^2 );
%
%  dV(1,2:n) = dV(1,2:n) - [ 1 ] * tau1 * ( [ 1; V(2:m,1) ]' * A(1:m,2:n) );       % 1.7
%  dV(1,2:n) = dV(1,2:n) - [ 1 ] * tau1 * ( A(1,2:n) + [ V(2:m,1) ]' * A(2:m,2:n) );
%  dV(1,2:n) = dV(1,2:n) - [ 1 ] * tau1 * ( A(1,2:n) + [ A(2:m,1)/ dalpha(1) ]' * A(2:m,2:n) );
%  dV(1,2:n) = dV(1,2:n) - [ 1 ] * tau1 * ( dV(1,2:n) + ( dG(1,2:n) + dA(2:n,1)' * dA(2:n,2:n) ) / dalpha(1) );
%  dV(1,2:n) = dV(1,2:n) - dtau(1) * ( dA(1,2:n) + ( dG(1,2:n) / dalpha(1) + dV(2:n,1)' * dA(2:n,2:n) )  );

   dV(1:n,2:n) = A(1:n,2:n) - [ [ 1; V(2:n,1) ] ] * dtau(1) * ( [ 1; V(2:m,1) ]' * A(1:m,2:n) );         % 1.7    --- isn't the point to not use V?

%
%  this is the check
   norm( triu(RR(1,1:n)) - triu(dV(1,1:n)), 'fro' ) / norm( triu(RR(1,1:n)), 'fro' )
   norm( tril(RR(1:n,1),-1) - tril(dV(1:n,1),-1), 'fro' ) / norm( tril(RR(1:n,1),-1), 'fro' )
   norm( V(1:n,1:n) - dV(1:n,1:n), 'fro' ) / norm( V(1:n,1:n), 'fro' ) 
%%%%%%
%
%  step 2
%  this is what we would do working on a V
   norma = norm(V(2:m,2),2);                                                                         % 2.1
   if ( V(2,2) > 0 ) V(2,2) = V(2,2) + norma; else V(2,2) = V(2,2) - norma; end                      % 2.2
   alpha(2) = V(2,2);                                                                                % 2.3
   V(3:m,2) = V(3:m,2) / alpha(2);                                                                   % 2.4
   if ( V(2,2) > 0 ) V(2,2) = - norma; else V(2,2) = norma; end                                      % 2.5
   tau2 = 2/( 1 + norm(V(3:m,2))^2 );                                                                % 2.6
   V(2:m,3:n) = A(2:m,3:n) - [ [ 1; V(3:m,2) ] ] * tau2 * ( [ 1; V(3:m,2) ]' * A(2:m,3:n) );         % 2.7
%   V(1:m,3:n) = A(1:m,3:n) - [ [0; 1; V(3:m,2) ] ] * tau2 * ( [0; 1; V(3:m,2) ]' * A(1:m,3:n) );         % 2.7
%
%  this is the check
   norm( triu(RR(2,1:n)) - triu(V(2,1:n)), 'fro' ) / norm( triu(RR(2,1:n)), 'fro' )
   norm( tril(RR(1:m,2),-1) - tril(V(1:m,2),-1), 'fro' ) / norm( tril(RR(1:m,2),-1), 'fro' )




return
%
%%%%%%
%
%  step 2
%  Update A'*A
%   
%   dG = RR(1,1:n)'*RR(1,1:n) - A'*A;                                              % 2.0 ??
%
%  norma = norm(V(1:m,1),2);                                                       % 2.1
   norma = sqrt( dA(1:n,2)'*dA(1:n,2) + dG(2,2) ); 
%
%  if ( V(1,1) > 0 ) V(1,1) = V(1,1) + norma; else V(1,1) = V(1,1) - norma; end    % 2.2
   if ( dV(2,2) > 0 ) dV(2,2) = dV(2,2) + norma; else dV(2,2) = dV(2,2) - norma; end
%
%  alpha(1) = V(1,1);                                                              % 2.3
   dalpha(2) = dV(2,2);
%
%  V(2:m,1) = V(2:m,1) / dalpha(1);                                                % 2.4
   dV(3:n,2) = dV(3:n,2) / dalpha(2);
%
%  if ( V(1,1) > 0 ) V(1,1) = - norma; else V(1,1) = norma; end                    % 2.5
   if ( dV(2,2) > 0 ) dV(2,2) = - norma; else dV(2,2) = norma; end
%
%  dtau(1) = 2/( 1 + norm(V(2:m,1))^2 );                                           % 2.6
   dtau(2) = 2/( 1 + dV(3:n,2)'*dV(3:n,2) + ( dG(2,2) ) / dalpha(2)^2 );
%
%  dV(1,2:n) = dV(1,2:n) - [ 1 ] * tau1 * ( [ 1; V(2:m,1) ]' * A(1:m,2:n) );       % 1.7
%  dV(1,2:n) = dV(1,2:n) - [ 1 ] * tau1 * ( A(1,2:n) + [ V(2:m,1) ]' * A(2:m,2:n) );
%  dV(1,2:n) = dV(1,2:n) - [ 1 ] * tau1 * ( A(1,2:n) + [ A(2:m,1)/ dalpha(1) ]' * A(2:m,2:n) );
%  dV(1,2:n) = dV(1,2:n) - [ 1 ] * tau1 * ( dV(1,2:n) + ( dG(1,2:n) + dA(2:n,1)' * dA(2:n,2:n) ) / dalpha(1) );
%  dV(1,2:n) = dV(1,2:n) - dtau(1) * ( dA(1,2:n) + ( dG(1,2:n) / dalpha(1) + dV(2:n,1)' * dA(2:n,2:n) )  );

   dV(2:n,3:n) = A(2:n,3:n) - [ [ 1; V(3:n,2) ] ] * dtau(2) * ( [ 1; V(3:m,2) ]' * A(2:m,3:n) );         % 2.7

%
%  this is the check
   norm( triu(RR(2,1:n)) - triu(dV(2,1:n)), 'fro' ) / norm( triu(RR(2,1:n)), 'fro' )
   norm( tril(RR(1:n,2),-1) - tril(dV(1:n,2),-1), 'fro' ) / norm( tril(RR(1:n,2),-1), 'fro' )
   norm( V(1:n,1:n) - dV(1:n,1:n), 'fro' ) / norm( V(1:n,1:n), 'fro' )

%%%%
