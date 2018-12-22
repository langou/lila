%
   function [ T ] = lila_larft_w0b_connect( n, i, j, A, T )
%
      T(i:j-1,j:j+n-1)   = ( - T(i:j-1,j:j+n-1) - T(i:j-1,1:j-1) * A(j:j+n-1,1:j-1)' ) ;
      T(i:j-1,j:j+n-1)   = T(i:j-1,j:j+n-1) / ( (eye(n,n) + tril(A(j:j+n-1,j:j+n-1),-1) )');
%
   end
