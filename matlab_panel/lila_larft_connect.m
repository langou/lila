%
   function [ T ] = lila_larft_connect( m, n, i, A, T )
%
      T(1:i-1,i:i+n-1) = ...
          - ( T(1:i-1,1:i-1) * A(i:m,1:i-1)' )* ...
              ( (tril(A(i:m,i:i+n-1),-1) + eye(m-i+1,n) ) * T(i:i+n-1,i:i+n-1) );
%
   end
