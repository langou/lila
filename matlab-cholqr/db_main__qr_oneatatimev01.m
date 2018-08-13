%
   fprintf('\n');
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
%%%%%%%
%
%  all quantites starting with a d must be computed with A only: dtau, dalpha, dbeta, dgamma, etc.
%
%%%%%%%
%
%  step 1   -  this is what we would do working on a V  
%
   norma1 = norm(V(1:m,1),2); 
   if ( V(1,1) > 0 ) V(1,1) = V(1,1) + norma1; else V(1,1) = V(1,1) - norma1; end
   V(2:m,1) = V(2:m) / V(1,1);
   if ( V(1,1) > 0 ) V(1,1) = - norma1; else V(1,1) = norma1; end
%
%  this is the check
%
   fprintf('Step 1 - Accuracy with V = %e', ...
   norm( triu(qr(As(1:m,1))) - [ V(1,1); zeros(m-1,1) ], 'fro' ) / norm( triu(qr(As(1:m,1))), 'fro' ) );
   fprintf('\n');
%
%  step 1   -  now we do the same (i.e. compute R(1,1)) without V  -- this one is easy because it is the first column
%
   dnorma1 = norm(A(1:m,1),2);
   if ( A(1,1) > 0 ) dalpha1 = A(1,1) + dnorma1; else dalpha1 = A(1,1) - dnorma1; end
   if ( A(1,1) > 0 ) dR(1,1) = - dnorma1; else dR(1,1) = dnorma1; end
%
%  this is the check 
%
   fprintf('Step 1 - Accuracy without V = %e', ...
   norm( triu(qr(As(1:m,1))) - [ dR(1,1); zeros(m-1,1) ], 'fro' ) / norm( triu(qr(As(1:m,1))), 'fro' ) );
   fprintf('\n');
%
%%%%%%%%%%
%
%  step 2     -   this is what we would do working on a V
%
   tau1 = 2/( 1+norm( V(2:m,1) )^2 );   tmp_product = [1;V(2:m,1)]' * A(1:m,2);
%
   V(1:m,2) = A(1:m,2) - [1;V(2:m,1)] * tau1 * tmp_product;
   V(2,2)   = A(2,2) -  V(2,1) * tau1 * tmp_product;
%
   norma2 = sqrt( V(2:m,2)'*V(2:m,2) );
%
   if ( V(2,2) > 0 ) V(2,2) = V(2,2) + norma2; else V(2,2) = V(2,2) - norma2; end
   V(3:m,2) = V(3:m,2) / V(2,2);
   if ( V(2,2) > 0 ) V(2,2) = - norma2; else V(2,2) = norma2; end
%
   fprintf('\n');
   fprintf('Step 2 - Accuracy with V = %e', ...
   norm( triu(qr(As(1:m,1:2))) - [ triu(V(1:2,1:2)); zeros(m-2,2) ], 'fro' ) / norm( triu(qr(As(1:m,1:2))), 'fro' ) );
   fprintf('\n');
%
%  step 2   -   now we do the same (i.e. compute R(1:2,1:2)) without V
%  step 2       this all the quantities from A that we will need:
%
%  we are only allowed to do a A'*A kind of stuff, so we do it here. 
%
   dbeta11 = A(2:m,1)' * A(2:m,1);                 
   dbeta12 = A(2:m,1)' * A(2:m,2);                
   dbeta22 = A(2:m,2)' * A(2:m,2);
%             
   dA12 = A(1,2);
   dA21 = A(2,1);
   dA22 = A(2,2);
%
%  The below, dividing by the previous alpha
%
   dtau1 = 2/( 1 + dbeta11 / dalpha1^2 ); 
%
   dgamma = A(1,2) + dbeta12 / dalpha1;              % this is [ 1; V(2:m,1) ]' * A(1:m,2)
%
   dR(1,2) = A(1,2) - dtau1 * dgamma;                % this is the correct R(1,2) computed only with A values :)
   dV(2,1) = A(2,1) / dalpha1;                       % I feel we will need this V(2,1) a lot (in step 3, step 4, etc.)
%
   dR(2,2) = A(2,2) - dV(2,1)* dtau1 * dgamma;       % this is not the final R
%
%%%%%%%
%  step 2: I am not sure how to compute dnorma2. Here is a bunch of ways to do it . . .
%
%  dnorma2 = norm( A(2:m,2) - [ A(2:m,1) ] * tau1 * dgamma / dalpha1 );                                            %  this might be a problem to compute
%
   dnorma2 = sqrt( dbeta22 - 2 * dbeta12 * dtau1 * dgamma / dalpha1  + dbeta11 *( dtau1 * dgamma / dalpha1 )^2 );  %  this might be a problem to compute accurately
%
%  dnorma2 = sqrt( norm(A(1:m,2))^2 - dR(1,2)^2 );                                                                 %  this is one way to do it, not sure if this is stable
%
%%%%%%
%
%  step 2: now, that we computed dnorma2 we can use it
%
   if ( dR(2,2) > 0 ) dR(2,2) = dR(2,2) + dnorma2; else dR(2,2) = dR(2,2) - dnorma2; end
   if ( dR(2,2) > 0 ) dR(2,2) = - dnorma2; else dR(2,2) = dnorma2; end
%
%  step 2: all done :) ready for the check
%
   fprintf('Step 2 - Accuracy without V = %e', ...
   norm( triu(qr(As(1:m,1:2))) - [ dR(1:2,1:2); zeros(m-2,2) ], 'fro' ) / norm( triu(qr(As(1:m,1:2))), 'fro' ) );
   fprintf('\n');

%%%%%%%%%%
return





