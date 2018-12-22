%
   clear
%
   m = 21;
   n = 6;
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
   [ RR ] = qr(A(1:m,1:n));
%
%  all quantites starting with a d must be computed with A only: dtau, dalpha, dbeta, dgamma, etc.
%
   dA(1:n,1:n) = A(1:n,1:n);
   dG(1:n,1:n) = A(n+1:m,1:n)' * A(n+1:m,1:n);
   dV(1:n,1:n) = dA(1:n,1:n);
   dGG(1:n,1:n) = A(1:m,1:n)' * A(1:m,1:n);
%
%===========================================================================================================================================
%
%  step 1
%  this is what we would do working on a V
   norma = norm(V(1:m,1),2);                                                                  % 1.1
   if ( V(1,1) > 0 ) V(1,1) = V(1,1) + norma; else V(1,1) = V(1,1) - norma; end               % 1.2
   alpha(1) = V(1,1);                                                                         % 1.3
   V(2:m,1) = V(2:m,1) / alpha(1);                                                            % 1.4
   if ( V(1,1) > 0 ) V(1,1) = - norma; else V(1,1) = norma; end                               % 1.5
   tau1 = 2/( 1 + norm(V(2:m,1))^2 );                                                         % 1.6
   V(1:m,2:n) = A(1:m,2:n) - [ [ 1; V(2:m,1) ] ] * tau1 * ( [ 1; V(2:m,1) ]' * A(1:m,2:n) );  % 1.7
%
%  this is the check
   fprintf('1 :: %e\n', norm( triu(RR(1,1:n)) - triu(V(1,1:n)), 'fro' ) / norm( triu(RR(1,1:n)), 'fro' ));
   fprintf('1 :: %e\n', norm( tril(RR(1:m,1),-1) - tril(V(1:m,1),-1), 'fro' ) / norm( tril(RR(1:m,1),-1), 'fro' ));
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
%  dV(1:n,2:n) = dV(1:n,2:n) - [ 1; V(2:n,1) ] * tau1 * ( [ 1; V(2:m,1) ]' * A(1:m,2:n) );       % 1.7
%  dV(1:n,2:n) = dV(1:n,2:n) - [ 1; V(2:n,1) ] * tau1 * ( A(1,2:n) + [ V(2:m,1) ]' * A(2:m,2:n) );
%  dV(1:n,2:n) = dV(1:n,2:n) - [ 1; dV(2:n,1) ] * tau1 * ( A(1,2:n) + [ A(2:m,1)/ dalpha(1) ]' * A(2:m,2:n) );
%  dV(1:n,2:n) = dV(1:n,2:n) - [ 1; dV(2:n,1) ] * tau1 * ( dV(1,2:n) + ( dG(1,2:n) + dA(2:n,1)' * dA(2:n,2:n) ) / dalpha(1) );
   dV(1:n,2:n) = dV(1:n,2:n) - [ 1; dV(2:n,1) ] * dtau(1) * ( dA(1,2:n) + ( dG(1,2:n) / dalpha(1) + dV(2:n,1)' * dA(2:n,2:n) )  );
%
%  this is the check
   fprintf('1 :: %e\n', norm( triu(RR(1,1:n)) - triu(dV(1,1:n)), 'fro' ) / norm( triu(RR(1,1:n)), 'fro' ));
   fprintf('1 :: %e\n', norm( tril(RR(1:n,1),-1) - tril(dV(1:n,1),-1), 'fro' ) / norm( tril(RR(1:n,1),-1), 'fro' ));
   fprintf('1 :: %e\n', norm( V(1:n,1:n) - dV(1:n,1:n), 'fro' ) / norm( V(1:n,1:n), 'fro' ));
%
%===========================================================================================================================================
%
%  step 2
%  this is what we would do working on a V
   V2(1:m,1:n) = V(1:m,1:n);                                                                         % save the current matrix so as to be able to rerun for check
   G(2:n,2:n) = V(n+1:m,2:n)'*V(n+1:m,2:n);                                                          % 2.0 ( compute the new normal equations from ``scratch``)
   GG(2:n,2:n) = V(2:m,2:n)'*V(2:m,2:n);                                                             % 2.0 ( compute the new normal equations from ``scratch``)
   norma = norm(V(2:m,2),2);                                                                         % 2.1
   if ( V(2,2) > 0 ) V(2,2) = V(2,2) + norma; else V(2,2) = V(2,2) - norma; end                      % 2.2
   alpha(2) = V(2,2);                                                                                % 2.3
   V(3:m,2) = V(3:m,2) / alpha(2);                                                                   % 2.4
   if ( V(2,2) > 0 ) V(2,2) = - norma; else V(2,2) = norma; end                                      % 2.5
   tau2 = 2/( 1 + norm(V(3:m,2))^2 );                                                                % 2.6
   V(2:m,3:n) = V(2:m,3:n) - [ [ 1; V(3:m,2) ] ] * tau2 * ( [ 1; V(3:m,2) ]' * V(2:m,3:n) );         % 2.7

%  this is the check
   fprintf('2.R :: %e\n',norm( triu(RR(2,1:n)) - triu(V(2,1:n)), 'fro' ) / norm( triu(RR(2,1:n)), 'fro' ));
   fprintf('2.R :: %e\n',norm( tril(RR(1:m,2),-1) - tril(V(1:m,2),-1), 'fro' ) / norm( tril(RR(1:m,2),-1), 'fro' ));
%
%  step 2.0
%   
   dG(2:n,2:n) = dG(2:n,2:n) + dA(1:n,2:n)' * dA(1:n,2:n) - dV(1:n,2:n)' * dV(1:n,2:n);             % 2.0
   dGG(1:n,1:n) = dGG(1:n,1:n) - dV(1,1:n)' * dV(1,1:n);                 
%
%  this is the check of step 2.0
%
%  norm ( ( dG(2:n,2:n) + dA(1:n,2:n)' * dA(1:n,2:n) - dV(1:n,2:n)' * dV(1:n,2:n) ) - ( V2(n+1:m,2:n)' * V2(n+1:m,2:n) ), 'fro' ) / norm( (V2(2:m,2:n)'*V2(2:m,2:n) ), 'fro' )
   fprintf('2.0 :: %e\n',norm ( dG(2:n,2:n) - G(2:n,2:n) , 'fro' ) / norm( G(2:n,2:n), 'fro' ));
   fprintf('2.0 :: %e\n',norm ( dGG(2:n,2:n) - GG(2:n,2:n) , 'fro' ) / norm( G(2:n,2:n), 'fro' ));
%
%  dnorma = norm( [ V2(2:m,2) ],2);                                                               % 2.1
   dnorma = sqrt( dGG(2,2) );                                                                     % 2.1
   fprintf('2.1 :: %e\n',abs( dnorma - norma ) / abs( norma ) );                                  % 2.1
%
%  if ( V(1,1) > 0 ) V(1,1) = V(1,1) + norma; else V(1,1) = V(1,1) - norma; end                   % 2.2
   dbeta(2) = dV(2,2);
   dGG(2,2) = dGG(2,2) - dV(2,2)^2;                                                               % 2.2
   if ( dV(2,2) > 0 ) dV(2,2) = dV(2,2) + norma; else dV(2,2) = dV(2,2) - norma; end
   dalpha(2) = dV(2,2);
   dGG(2,2) = dGG(2,2) + dV(2,2)^2;
%
%  V(3:m,2) = V(3:m,2) / alpha(2);                                                                % 2.4
   dV(3:n,2) = dV(3:n,2) / dalpha(2);                                                             % 2.4
%
%  if ( V(2,2) > 0 ) V(2,2) = - norma; else V(2,2) = norma; end                                   % 2.5
   if ( dV(2,2) > 0 ) dV(2,2) = - norma; else dV(2,2) = norma; end                                % 2.5
%
%  dtau(2) = 2/( 1 + norm(V(3:m,2))^2 );                                                          % 2.6
   dtau(2) = 2 / ( ( dGG(2,2) ) / dalpha(2)^2 );                                                  % 2.6
   fprintf('2.6 :: %e\n',abs( dtau(2) - tau2 ) / abs( tau2 ) );                                   % 2.6
%
%  dV(2:n,3:n) = V2(2:n,3:n) - [ 1; V(3:n,2) ] * dtau(2) * ( [ 1; V(3:m,2) ]' * V2(2:m,3:n) );                                   % 2.7
%  dV(2:n,3:n) = dV(2:n,3:n) - [ 1; dV(3:n,2) ] * dtau(2) * ( [ 1; V(3:m,2) ]' * V2(2:m,3:n) );                                  % 2.7
%  dV(2:n,3:n) = dV(2:n,3:n) - [ 1; dV(3:n,2) ] * dtau(2) * ( V2(2,3:n) + V(3:m,2)' * V2(3:m,3:n) );                             % 2.7
%  dV(2:n,3:n) = dV(2:n,3:n) - [ 1; dV(3:n,2) ] * dtau(2) * ( V2(2,3:n) + ( dGG(2,3:n) - V2(2,2)' * V2(2,3:n) ) / dalpha(2) );   % 2.7
   dV(2:n,3:n) = dV(2:n,3:n) - [ 1; dV(3:n,2) ] * dtau(2) * ( dV(2,3:n) + ( dGG(2,3:n) - dbeta(2)' * dV(2,3:n) ) / dalpha(2) );  % 2.7
%
%  this is the check
   fprintf('2.7 :: %e\n',norm( triu(RR(2,1:n)) - triu(dV(2,1:n)), 'fro' ) / norm( triu(RR(2,1:n)), 'fro' ));
   fprintf('2.7 :: %e\n',norm( tril(RR(1:n,2),-1) - tril(dV(1:n,2),-1), 'fro' ) / norm( tril(RR(1:n,2),-1), 'fro' ));
   fprintf('2.7 :: %e\n',norm( V(1:n,1:n) - dV(1:n,1:n), 'fro' ) / norm( V(1:n,1:n), 'fro' ));


return

%
%===========================================================================================================================================
%
   k = 3;
%
%===========================================================================================================================================
%
%  step 3
%  this is what we would do working on a V
   V3(1:m,1:n) = V(1:m,1:n);                                                                                   % save the current matrix so as to be able to rerun for check
   G1(k:n,k:n) = V(n+1:m,k:n)'*V(n+1:m,k:n);                                                                   % 3.0 ( compute the new normal equations from ``scratch``)
   GG1(k:n,k:n) = V(k:m,k:n)'*V(k:m,k:n); 
   norma = norm(V(k:m,k),2);                                                                                   % 3.1
   if ( V(k,k) > 0 ) V(k,k) = V(k,k) + norma; else V(k,k) = V(k,k) - norma; end                                % 3.2
   alpha(k) = V(k,k);                                                                                          % 3.3
   V(k+1:m,k) = V(k+1:m,k) / alpha(k);                                                                         % 3.4
   if ( V(k,k) > 0 ) V(k,k) = - norma; else V(k,k) = norma; end                                                % 3.5 
   tau3 = 2/( 1 + norm(V(k+1:m,k))^2 );                                                                        % 3.6
   V(k:m,k+1:n) = V(k:m,k+1:n) - [ [ 1; V(k+1:m,k) ] ] * tau3 * ( [ 1; V(k+1:m,k) ]' * V(k:m,k+1:n) );         % 3.7

%  this is the check
   fprintf('1 :: %e\n',norm( triu(RR(k,1:n)) - triu(V(k,1:n)), 'fro' ) / norm( triu(RR(k,1:n)), 'fro' ));
   fprintf('1 :: %e\n',norm( tril(RR(1:m,k),-1) - tril(V(1:m,k),-1), 'fro' ) / norm( tril(RR(1:m,k),-1), 'fro' ));
%
%%%%
%
%  step 3.0
%  
   dGG(2:n,2:n) = dGG(2:n,2:n) - V(2,2:n)' * V(2,2:n);                 



return


    

    dG(3:n,3:n) = G1(3:n,3:n);   %This works so I am calculating the wrong normal equation
    dGG(3:n,3:n) = G1(3:n,3:n);  %This works so I am calculating the wrong normal equation



%  this is the check of step 2.0
%
   norm ( dG(k:n,k:n) - G1(k:n,k:n) , 'fro' ) / norm( G1(k:n,k:n), 'fro' )
   norm ( dGG(k:n,k:n) - G1(k:n,k:n) , 'fro' ) / norm( G1(k:n,k:n), 'fro' )

%
%  step 2.1
%
%  dnorma = norm( [ V2(2:m,2) ],2);                                                                          % 3.1
%  dnorma = sqrt( [ V2(2:m,2)]'*[V2(2:m,2) ] );                                                              % 3.1
%  dnorma = sqrt( V2(2:n,2)'*V2(2:n,2) + [ V2(n+1:m,2)]'*[V2(n+1:m,2) ] );                                   % 3.1
%  dnorma = sqrt( dV(2:n,2)'*dV(2:n,2) + [ V2(n+1:m,2)]'*[V2(n+1:m,2) ] );                                   % 3.1
%   dnorma = sqrt( dV(2:n,2)'*dV(2:n,2) + dG(2,2) );                                                         % 3.1
%   dnorma = sqrt( dV(k:n,k)'*dV(k:n,k) + dG(k,k) );                                                          % 3.1
   dnorma = sqrt( dV(k:n,k)'*dV(k:n,k) + dGG(k,k) );                                                          % 3.1
%
%  this is the check of step 2.1
%
   abs( dnorma - norma ) / abs( norma )
%
%  if ( V(1,1) > 0 ) V(1,1) = V(1,1) + norma; else V(1,1) = V(1,1) - norma; end                             % 3.2
%   if ( dV(2,2) > 0 ) dV(2,2) = dV(2,2) + norma; else dV(2,2) = dV(2,2) - norma; end
   if ( dV(k,k) > 0 ) dV(k,k) = dV(k,k) + norma; else dV(k,k) = dV(k,k) - norma; end
%
%  alpha(1) = V(1,1);                                                                                       % 3.3
   dalpha(2) = dV(2,2);
   dalpha(k) = dV(k,k);
%   dalpha(2) = alpha(2);
%
%  V(2:m,1) = V(2:m,1) / dalpha(1);                                                                         % 3.4
%   dV(3:n,2) = dV(3:n,2) / dalpha(2);
   dV(k+1:n,k) = dV(k+1:n,k) / dalpha(k);
%
%  if ( V(1,1) > 0 ) V(1,1) = - norma; else V(1,1) = norma; end                                             % 3.5
%   if ( dV(2,2) > 0 ) dV(2,2) = - norma; else dV(2,2) = norma; end
   if ( dV(k,k) > 0 ) dV(k,k) = - norma; else dV(k,k) = norma; end
%
%  dtau(1) = 2/( 1 + norm(V(2:m,1))^2 );                                                                    % 3.6
%   dtau(2) = 2/( 1 + dV(3:n,2)'*dV(3:n,2) + ( dG(2,2) ) / dalpha(2)^2 );
%   dtau(k) = 2/( 1 + dV(k+1:n,k)'*dV(k+1:n,k) + ( dG(k,k) ) / dalpha(k)^2 );
   dtau(k) = 2/( 1 + dV(k+1:n,k)'*dV(k+1:n,k) + ( dGG(k,k) ) / dalpha(k)^2 );
%
%  dV(2:n,3:n) = dV(2:n,3:n) - [ 1; V(3:n,2) ] * tau2 * ( [ 1; V(3:m,2) ]' * A(2:m,3:n) );       % 1.7
%  dV(2:n,3:n) = dV(2:n,3:n) - [ 1; V(3:n,2) ] * tau2 * ( A(2,3:n) + [ V(3:m,2) ]' * A(3:m,3:n) );
%  dV(2:n,3:n) = dV(2:n,3:n) - [ 1; dV(3:n,2) ] * tau2 * ( A(2,3:n) + [ A(3:m,2)/ dalpha(3) ]' * A(3:m,3:n) );
%  dV(2:n,3:n) = dV(2:n,3:n) - [ 1; dV(3:n,2) ] * tau2 * ( dV(2,3:n) + ( dG(2,3:n) + dA(3:n,2)' * dA(3:n,3:n) ) / dalpha(2) );
%  dV(2:n,3:n) = dV(2:n,3:n) - [ 1; dV(3:n,2) ] * dtau(2) * ( dA(2,3:n) + ( dG(2,3:n) / dalpha(2) + dV(3:n,2)' * dA(3:n,3:n) )  );
%
%   dV(2:n,3:n) = dV(2:n,3:n) - [ 1; dV(3:n,2) ] * dtau(2) * ( dV(2,3:n) + ( dG(2,3:n) / dalpha(2) + dV(3:n,2)' * dV(3:n,3:n) )  );      % 3.7
%   dV(k:n,k+1:n) = dV(k:n,k+1:n) - [ 1; dV(k+1:n,k) ] * dtau(k) * ( dV(k,k+1:n) + ( dG(k,k+1:n) / dalpha(k) + dV(k+1:n,k)' * dV(k+1:n,k+1:n) )  );      % 3.7
   dV(k:n,k+1:n) = dV(k:n,k+1:n) - [ 1; dV(k+1:n,k) ] * dtau(k) * ( dV(k,k+1:n) + ( dGG(k,k+1:n) / dalpha(k) + dV(k+1:n,k)' * dV(k+1:n,k+1:n) )  );      % 3.7
%
%  this is the check
%   norm( triu(RR(2,1:n)) - triu(dV(2,1:n)), 'fro' ) / norm( triu(RR(2,1:n)), 'fro' )
   fprintf('1 :: %e\n',norm( triu(RR(k,1:n)) - triu(dV(k,1:n)), 'fro' ) / norm( triu(RR(k,1:n)), 'fro' ));
%   norm( tril(RR(1:n,2),-1) - tril(dV(1:n,2),-1), 'fro' ) / norm( tril(RR(1:n,2),-1), 'fro' )
   fprintf('1 :: %e\n',norm( tril(RR(1:n,k),-1) - tril(dV(1:n,k),-1), 'fro' ) / norm( tril(RR(1:n,k),-1), 'fro' ));
%   norm( V(1:n,1:n) - dV(1:n,1:n), 'fro' ) / norm( V(1:n,1:n), 'fro' )
   fprintf('1 :: %e\n',norm( V(1:n,1:n) - dV(1:n,1:n), 'fro' ) / norm( V(1:n,1:n), 'fro' ));
%
%===========================================================================================================================================
%
