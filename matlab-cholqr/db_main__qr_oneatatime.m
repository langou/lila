%
   clear
%
   m = 20;
   n = 4;
   log10KA = 3;
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
%  all quantites starting with a d must be computed with A only: dtau, dalpha, dbeta, dgamma, etc.
%
%  step 1
%  this is what we would do working on a V
   norma1 = norm(V(1:m,1),2); 
   if ( V(1,1) > 0 ) V(1,1) = V(1,1) + norma1; else V(1,1) = V(1,1) - norma1; end
   alpha1 = V(1,1);
   V(2:m,1) = V(2:m) / V(1,1);
   if ( V(1,1) > 0 ) V(1,1) = - norma1; else V(1,1) = norma1; end
%  this is the check
   norm( triu(qr(As(1:m,1))) - [ V(1,1); zeros(m-1,1) ], 'fro' ) / norm( triu(qr(As(1:m,1))), 'fro' )
%
%  step 1
%  now we do the same (i.e. compute R(1,1)) without V
   dnorma1 = norm(A(1:m,1),2);
   if ( A(1,1) > 0 ) dalpha1 = A(1,1) + dnorma1; else dalpha1 = A(1,1) - dnorma1; end
   if ( A(1,1) > 0 ) dR(1,1) = - dnorma1; else dR(1,1) = dnorma1; end
%  this is the check
   norm( triu(qr(As(1:m,1))) - [ dR(1,1); zeros(m-1,1) ], 'fro' ) / norm( triu(qr(As(1:m,1))), 'fro' )
%  now we do the same (i.e. compute R(1,1)) without V
%
%  step 2
%  this is what we would do working on a V
   tau1 = 2/( 1 + norm(V(2:m,1))^2 );
   V(1:m,2) = A(1:m,2) - [ 1; V(2:m,1) ] * tau1 * ( [ 1; V(2:m,1) ]' * A(1:m,2) );
   norma2 = norm(V(2:m,2),2);
   if ( V(2,2) > 0 ) V(2,2) = V(2,2) + norma2; else V(2,2) = V(2,2) - norma2; end
   V(3:m,2) = V(3:m,2) / V(2,2);
   if ( V(2,2) > 0 ) V(2,2) = - norma2; else V(2,2) = norma2; end
%  this is the check
   norm( triu(qr(As(1:m,1:2))) - [ triu(V(1:2,1:2)); zeros(m-2,2) ], 'fro' ) / norm( triu(qr(As(1:m,1:2))), 'fro' )
%
%  step 2
%  now we do the same (i.e. compute R(1:2,1:2)) without V
   dtau1 = 2/( 1 + norm(A(2:m,1))^2/dalpha1^2 ); 
   dbeta11 = A(2:m,1)' * A(2:m,1);        % we are only allowed to do a A'*A kind of stuff, so we do it here. 
   dbeta12 = A(2:m,1)' * A(2:m,2);        % we are only allowed to do a A'*A kind of stuff, so we do it here. 
   dbeta22 = A(2:m,2)'* A(2:m,2);
   dgamma = A(1,2) + dbeta12 / dalpha1;   % this is [ 1; V(2:m,1) ]' * A(1:m,2)
   dR(1,2) = A(1,2) - dtau1 * dgamma;     % this is the correct R(1,2) computed only with A values :)
   dV(2,1) = A(2,1) / dalpha1;
   dR(2,2) = A(2,2) - [ dV(2,1) ] * dtau1 * dgamma;
%  I am not sure how to compute dnorma2. Here is a bunch of way to do it . . .
%  dnorma2 = norm( A(2:m,2) - [ A(2:m,1) ] * tau1 * dgamma / dalpha1 );                                            %  this might be a problem to compute
   dnorma2 = sqrt( dbeta22 - 2 * dbeta12 * tau1 * dgamma / dalpha1  + dbeta11 *( tau1 * dgamma / dalpha1 )^2 );    %  this might be a problem to compute accurately
%  dnorma2 = sqrt( norm(A(1:m,2))^2 - dR(1,2)^2 );                                                                 %  this is one way to do it, not sure if this is stable
   if ( dR(2,2) > 0 ) dR(2,2) = dR(2,2) + dnorma2; else dR(2,2) = dR(2,2) - dnorma2; end
   if ( dR(2,2) > 0 ) dR(2,2) = - dnorma2; else dR(2,2) = dnorma2; end
   norm( triu(qr(As(1:m,1:2))) - [ dR(1:2,1:2); zeros(m-2,2) ], 'fro' ) / norm( triu(qr(As(1:m,1:2))), 'fro' )





return









%
%
%
%  alpha2 = norm(A(2:m,2),2);
%  a2 = A(1:m,2) - a1 *( tau1*(a1'*A(1:m,2)) );
%  a2(1) = a2(1) + alpha2;
%  a2 = a2/alpha2;
%  tau2 = 2/(1+norm(a2,2))^2;
%
%
%
%  alpha3 = norm(A(1:m,3),2);
%  a3 = A(1:m,3) - a2 * ( tau2*(a2'*A(1:m,3)) ) - a1 * ( tau1*(a1'*A(1:m,3)) ) + a2 * ( tau2*(a2'*a1)*tau1*(a1'*A(1:m,3)) );
%  a3(1) = a3(1) + alpha3;
%  a3 = a3/alpha3;
%  tau3 = 2/(1+norm(a3,2))^2;
%
%
%
%  alpha4 = norm(A(1:m,4),2);
%  a4 = A(1:m,4) - a2 * ( tau2*(a2'*A(1:m,4)) ) - a1 * ( tau1*(a1'*A(1:m,4)) ) + a2 * ( tau2*(a2'*a1)*tau1*(a1'*A(1:m,4)) ) - a3 * (tau3*(a3'*A(1:m,4))) + a3 * (tau3*(a3'*a2)*tau2*(a2'*A(1:m,4))) + a3 * (tau3*(a3'*a1)*tau1*(a1'*A(1:m,4))) - a3 * (tau3*(a3'*a2)*(tau2)*(a2'*a1)*(tau1)*(a1'*A(1:m,4)) );
%  a4(1) = a4(1) + alpha4;
%  a4 = a4/alpha4;
%  tau4 = 2/(1+norm(a4,2))^2;
%
   
