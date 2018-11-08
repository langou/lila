%
   function [ A, T, Q ] = lila_geqr2_w0b( m, n, A, T, Q )
%
      Q(1:m,1:n) = A(1:m,1:n);
%
%     Cholesky QR
%
      R(1:n,1:n) = Q(1:m,1:n)'*Q(1:m,1:n); 
      R(1:n,1:n) = chol( R(1:n,1:n), 'upper' ); 
      Q(1:m,1:n) = Q(1:m,1:n) / R(1:n,1:n);
%
%     Householder reconstruction
%
      A(1:m,1:n) = Q(1:m,1:n); 
      D = -eye(n,n);
%
      for k = 1:n,
         if (abs(1 - A(k,k)) < abs( 1 + A(k,k) ))
            Q(1:m,k) = -Q(1:m,k);
            R(k,k:n) = -R(k,k:n);
         else 
            A(1:m,k) = -A(1:m,k);
            D(k,k) = 1.0e+00;
         end
         A(k,k) = 1 + A(k,k);
         A(1+k:m,k) = A(1+k:m,k) / A(k,k);
         A(1+k:m,1+k:n) = A(1+k:m,1+k:n) - A(1+k:m,k) * A(k,1+k:n);
      end
%
%    Construction of T
%
      T(1:n,1:n) = triu( A(1:n,1:n) ) / ( (eye(n,n) + tril(A(1:n,1:n),-1) )');
%
%    Put back the good R
%
      for ii=1:n,
         for jj=ii:n,
           A(ii,jj) = R(ii,jj);
         end
      end
%
   end

