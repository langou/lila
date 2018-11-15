%
   function [ A, T, Q, D ] = lila_orgh2_w0b_panel( m, n, j, A, T, Q, D )
%
      for k = 1:n,
         if (abs( 1 - T(j+k-1,j+k-1) ) < abs( 1 + T(j+k-1,j+k-1) ))
           T(j:j+n-1,j+k-1)     = - T(j:j+n-1,j+k-1);
           D(j+k-1)             = -1;
         else
           D(j+k-1)             = 1;
         end
         T(j+k-1,j+k-1)         = T(j+k-1,j+k-1) - 1;
         T(j+k:j+n-1,j+k-1)     = T(j+k:j+n-1,j+k-1) / T(j+k-1,j+k-1);
         T(j+k:j+n-1,j+k:j+n-1) = T(j+k:j+n-1,j+k:j+n-1) - T(j+k:j+n-1,j+k-1) * T(j+k-1,j+k:j+n-1);
      end
%
      for k = 1:n,
         if (D(j+k-1) == -1),
           A(j+n:m,j+k-1) = - A(j+n:m,j+k-1);
         end
      end
%
      A(j+n:m,j:j+n-1) = A(j+n:m,j:j+n-1) / triu(T(j:j+n-1,j:j+n-1));
%
      for k = 1:n,
         if (D(j+k-1) == -1),
           T(1:j-1,j+k-1) = - T(1:j-1,j+k-1);
         end
      end
%
      for k = 1:n,
         if (D(j+k-1) == -1),
           A(j+k-1,j+k-1:j+n-1) = - A(j+k-1,j+k-1:j+n-1);
         end
      end
%
m
      for k = 1:n,
         if (D(j+k-1) == -1),
           Q(1:m,j+k-1) = - Q(1:m,j+k-1);
         end
      end
%
      for ii=1:n,
         for jj=1:ii-1,
           A(j+ii-1,j+jj-1) = T(j+ii-1,j+jj-1);
         end
      end
%
      T(j:j+n-1,j:j+n-1) = - triu( T(j:j+n-1,j:j+n-1) ) / ( (eye(n,n) + tril(A(j:j+n-1,j:j+n-1),-1) )');
%
   end