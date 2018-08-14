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
   dA(1,1) = A(1,1);
   dC(1,1) = A(2:m,1)' * A(2:m,1);                 % we are only allowed to do a A'*A kind of stuff, so we do it here. 
%
   dnorma = sqrt( dC(1,1) + dA(1,1)*dA(1,1) );
   dR(1,1) = dA(1,1);
   if ( dR(1,1) > 0 ) dR(1,1) = dR(1,1) + dnorma; else dR(1,1) = dR(1,1) - dnorma; end
   dalpha(1) = dR(1,1);
   if ( dR(1,1) > 0 ) dR(1,1) = - dnorma; else dR(1,1) = dnorma; end
%  this is the check
   norm( triu(qr(As(1:m,1))) - [ dR(1,1); zeros(m-1,1) ], 'fro' ) / norm( triu(qr(As(1:m,1))), 'fro' )
%  now we do the same (i.e. compute R(1,1)) without V
%
%  step 2
%  this is what we would do working on a V
   tau1 = 2/( 1 + norm(V(2:m,1))^2 );
   V(1:m,2) = A(1:m,2) - [ 1; V(2:m,1) ] * tau1 * ( [ 1; V(2:m,1) ]' * A(1:m,2) );
   norma2_3 = V(3:m,2)'*V(3:m,2);
   norma2 = sqrt( V(2,2)*V(2,2) + norma2_3 );
   if ( V(2,2) > 0 ) V(2,2) = V(2,2) + norma2; else V(2,2) = V(2,2) - norma2; end
   V(3:m,2) = V(3:m,2) / V(2,2);
   if ( V(2,2) > 0 ) V(2,2) = - norma2; else V(2,2) = norma2; end
%  this is the check
   norm( triu(qr(As(1:m,1:2))) - [ triu(V(1:2,1:2)); zeros(m-2,2) ], 'fro' ) / norm( triu(qr(As(1:m,1:2))), 'fro' )
%
%  step 2
%  now we do the same (i.e. compute R(1:2,1:2)) without V
%
%  step 2: this all the quantities from A that we will need:
   dbeta11_2 = A(2:m,1)' * A(2:m,1);                 % we are only allowed to do a A'*A kind of stuff, so we do it here. 
   dbeta12_2 = A(2:m,1)' * A(2:m,2);                 % we are only allowed to do a A'*A kind of stuff, so we do it here. 
   dbeta22_2 = A(2:m,2)' * A(2:m,2);                 % we are only allowed to do a A'*A kind of stuff, so we do it here.
   dA12 = A(1,2);
   dA21 = A(2,1);
   dA22 = A(2,2);
%
   dtau1 = 2/( 1 + dbeta11_2 / dalpha(1)^2 ); 
   dgamma = dA12 + dbeta12_2 / dalpha(1);            % this is [ 1; V(2:m,1) ]' * A(1:m,2)
   dR(1,2) = dA12 - dtau1 * dgamma;                % this is the correct R(1,2) computed only with A values :)
   dV(2,1) = dA21 / dalpha(1);                     % I feel we will need this V(2,1) a lot (in step 3, step 4, etc.)
   dR(2,2) = dA22 - [ dV(2,1) ] * dtau1 * dgamma;  % this is not the final R
   dR22 = dR(2,2);
%
%  step 2: I am not sure how to compute dnorma2. Here is a bunch of way to do it . . .
%  dnorma2 = norm( A(2:m,2) - [ A(2:m,1) ] * tau1 * dgamma / dalpha(1) );                                              %  this might be a problem to compute
   dnorma2 = sqrt( dbeta22_2 - 2 * dbeta12_2 * dtau1 * dgamma / dalpha(1)  + dbeta11_2 *( dtau1 * dgamma / dalpha(1) )^2 );  %  this might be a problem to compute accurately
%  dnorma2 = sqrt( norm(A(1:m,2))^2 - dR(1,2)^2 );                                                                     %  this is one way to do it, not sure if this is stable
%
%  step 2: now, that we computed dnorma2 we can use it
   if ( dR(2,2) > 0 ) dR(2,2) = dR(2,2) + dnorma2; else dR(2,2) = dR(2,2) - dnorma2; end
   dalpha(2) = dR(2,2);
   if ( dR(2,2) > 0 ) dR(2,2) = - dnorma2; else dR(2,2) = dnorma2; end
%
%  step 2: all done :) ready for the check
   norm( triu(qr(As(1:m,1:2))) - [ dR(1:2,1:2); zeros(m-2,2) ], 'fro' ) / norm( triu(qr(As(1:m,1:2))), 'fro' )
%
%  step 3
%  this is what we would do working on a V
   tau2 = 2/( 1 + norm(V(3:m,2))^2 );
   V(1:m,3) = A(1:m,3);
   V(1:m,3) = V(1:m,3) - [ 1; V(2:m,1) ] * tau1 * ( [ 1; V(2:m,1) ]' * V(1:m,3) );
   V(1:m,3) = V(1:m,3) - [ 0; 1; V(3:m,2) ] * tau2 * ( [ 0; 1; V(3:m,2) ]' * V(1:m,3) );
   norma3 = norm(V(3:m,3),2);
   if ( V(3,3) > 0 ) V(3,3) = V(3,3) + norma3; else V(3,3) = V(3,3) - norma3; end
   V(4:m,3) = V(4:m,3) / V(3,3);
   if ( V(3,3) > 0 ) V(3,3) = - norma3; else V(3,3) = norma3; end
%  this is the check
   norm( triu(qr(As(1:m,1:3))) - [ triu(V(1:3,1:3)); zeros(m-3,3) ], 'fro' ) / norm( triu(qr(As(1:m,1:3))), 'fro' )
%
%  step 2
%  now we do the same (i.e. compute R(1:2,1:2)) without V
%
%  step 2: this all the quantities from A that we will need:
%
   dbeta22_3 = A(3:m,2)' * A(3:m,2);                 % we are only allowed to do a A'*A kind of stuff, so we do it here.
   dbeta12_3 = A(3:m,1)' * A(3:m,2);                 % we are only allowed to do a A'*A kind of stuff, so we do it here.
   dbeta11_3 = A(3:m,1)' * A(3:m,1);                 % we are only allowed to do a A'*A kind of stuff, so we do it here.
   dbeta12_2 = A(2:m,1)' * A(2:m,2);                 % we are only allowed to do a A'*A kind of stuff, so we do it here.
%
%  VV = norm( A(3:m,2) - [ V(3:m,1) ] * tau1 * ( [ 1; V(2:m,1) ]' * A(1:m,2) ) );
   VV = dbeta22_3 - 2* dtau1 * ( dA12 + [ dbeta12_2/dalpha(1) ] ) * ( dbeta12_3/dalpha(1) )  + dbeta11_3/dalpha(1)^2 * ( dtau1 * ( dA12 + [ dbeta12_2/dalpha(1) ] ) )^2  ;
   dtau2 = 2/( 1 + VV/( dalpha(2)^2 ) );
   abs(dtau2 - tau2 )
%
%
%
   norma3 = norm(V(3:m,3),2);
   VV(1:m,3) = A(1:m,3);
   VV(2:m,3) = VV(2:m,3) - [ V(2:m,1) ] * dtau1 * ( [ 1; V(2:m,1) ]' * VV(1:m,3) );
   VV(3:m,3) = VV(3:m,3) - [ V(3:m,2) ] * dtau2 * ( [ 1; V(3:m,2) ]' * VV(2:m,3) );
   dnorma3 = norm(V(3:m,3),2);
   abs(dnorma3 - norma3 )
%
   norma3 = norm(V(3:m,3),2);
   VV(2:m,3) = VV(2:m,3) - [ V(2:m,1) ] * dtau1 * ( [ 1; V(2:m,1) ]' * A(1:m,3) );
   VV(3:m,3) = VV(3:m,3) - [ V(3:m,2) ] * dtau2 * ( [ 1; V(3:m,2) ]' * VV(2:m,3) );
   dnorma3 = norm(V(3:m,3),2);
   abs(dnorma3 - norma3 )
%












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
   
