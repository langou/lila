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
%
%
%
   alpha1 = norm(A(1:m,1),2); 
   a1 = A(1:m,1);
   a1(1) = a1(1) + alpha1;
   a1 = a1/alpha1;
   tau1 = 2/(1+norm(a1,2))^2;
%
%
%
   alpha2 = norm(A(2:m,2),2);
   a2 = A(1:m,2) - a1 *( tau1*(a1'*A(1:m,2)) );
   a2(1) = a2(1) + alpha2;
   a2 = a2/alpha2;
   tau2 = 2/(1+norm(a2,2))^2;
%
%
%
   alpha3 = norm(A(1:m,3),2);
   a3 = A(1:m,3) - a2 * ( tau2*(a2'*A(1:m,3)) ) - a1 * ( tau1*(a1'*A(1:m,3)) ) + a2 * ( tau2*(a2'*a1)*tau1*(a1'*A(1:m,3)) );
   a3(1) = a3(1) + alpha3;
   a3 = a3/alpha3;
   tau3 = 2/(1+norm(a3,2))^2;
%
%
%
   alpha4 = norm(A(1:m,4),2);
   a4 = A(1:m,4) - a2 * ( tau2*(a2'*A(1:m,4)) ) - a1 * ( tau1*(a1'*A(1:m,4)) ) + a2 * ( tau2*(a2'*a1)*tau1*(a1'*A(1:m,4)) ) - a3 * (tau3*(a3'*A(1:m,4))) + a3 * (tau3*(a3'*a2)*tau2*(a2'*A(1:m,4))) + a3 * (tau3*(a3'*a1)*tau1*(a1'*A(1:m,4))) - a3 * (tau3*(a3'*a2)*(tau2)*(a2'*a1)*(tau1)*(a1'*A(1:m,4)) );
   a4(1) = a4(1) + alpha4;
   a4 = a4/alpha4;
   tau4 = 2/(1+norm(a4,2))^2;
%
%
%
   for j=1:n,

      A1 = As(1:m,1:j-1)'*As(1:m,1:j-1);
      T1 = inv( A1 );
      T2 = 2*T1 - T1*A1*T1;

      w = (eye(m) - As(1:m,1:j-1) * T2 * As(1:m,1:j-1)' ) * As(1:m,j);
      Q(:,j) = w / norm(w);

   end
   R = Q'*As;
%
   orth = norm(eye(n) - Q'*Q, 'fro');
   repres = norm( As - Q * R, 'fro') / norm( As, 'fro' ) ;
   fprintf('|| I - Q''* Q ||           = %6.1e\n', orth );
   fprintf('|| A - Q * R || / || A || = %6.1e\n', repres );
%

