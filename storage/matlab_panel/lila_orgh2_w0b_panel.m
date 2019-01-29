%
   function [ A, T, Q, D ] = lila_orgh2_w0b_panel( m, n, i, j, A, T, Q, D )
%
      for k = 1:n,
%
         if (abs(T(j+k-1,j+k-1) - 1) < abs( - T(j+k-1,j+k-1) - 1 ))
           Q(i:m,j+k-1)         = - Q(i:m,j+k-1);
           A(j+k-1,j+k-1:j+n-1) = - A(j+k-1,j+k-1:j+n-1);
           A(j+n:m,j+k-1)       = - A(j+n:m,j+k-1);
           T(j:j+n-1,j+k-1)     = - T(j:j+n-1,j+k-1);
           T(i:j-1,j+k-1)       = - T(i:j-1,j+k-1);
           D(j+k-1)             = -1;
         else
           D(j+k-1)             = 1;
         end
         T(j+k-1,j+k-1)         = T(j+k-1,j+k-1) - 1;
         T(j+k:j+n-1,j+k-1)     = T(j+k:j+n-1,j+k-1) / T(j+k-1,j+k-1);
         A(j+n:m,j+k-1)         = A(j+n:m,j+k-1) / T(j+k-1,j+k-1);
         T(j+k:j+n-1,j+k:j+n-1) = T(j+k:j+n-1,j+k:j+n-1) - T(j+k:j+n-1,j+k-1) * T(j+k-1,j+k:j+n-1);
         A(j+n:m,j+k:j+n-1)     = A(j+n:m,j+k:j+n-1) - A(j+n:m,j+k-1) * T(j+k-1,j+k:j+n-1);
      end
%
      for ii=1:n,
         for jj=1:ii-1,
           A(j+ii-1,j+jj-1) = T(j+ii-1,j+jj-1);
         end
      end
%
      T(j:j+n-1,j:j+n-1) = - triu( T(j:j+n-1,j:j+n-1) ) / ( (eye(n,n) + tril(A(j:j+n-1,j:j+n-1),-1) )');

   end

