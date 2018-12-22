%
   function [ Q, R ] = Cholesky_qr( Q )
%
      [ m, n ] = size(Q);
%
      R(1:n,1:n) = Q(1:m,1:n)'*Q(1:m,1:n); 
      R(1:n,1:n) = chol( R(1:n,1:n), 'upper' ); 
      Q(1:m,1:n) = Q(1:m,1:n) / R(1:n,1:n);
%
   end

