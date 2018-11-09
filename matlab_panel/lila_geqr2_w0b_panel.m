%
   function [ A, T, Q, D ] = lila_geqr2_w0b_panel( m, n, i, A, T, Q, D )
%
      for k = 1:i-1,
         if (D(k) == -1)
           A(k,i:i+n-1) = - A(k,i:i+n-1);
         end
      end
%
      T(1:i-1,i:i+n-1)   = Q(1:i-1,i:i+n-1);
      T(i:i+n-1,i:i+n-1) = Q(i:i+n-1,i:i+n-1);
      A(i+n:m,i:i+n-1)   = Q(i+n:m,i:i+n-1); 
%
      if( i ~= 1),
         T(1:i-1,i:i+n-1) =  (tril(A(1:i-1,1:i-1),-1) + eye(i-1,i-1)) \ T(1:i-1,i:i+n-1);
         T(i:i+n-1,i:i+n-1)  = T(i:i+n-1,i:i+n-1) - A(i:i+n-1,1:i-1) * T(1:i-1,i:i+n-1);
         A(i+n:m,i:i+n-1)  = A(i+n:m,i:i+n-1) - A(i+n:m,1:i-1) * T(1:i-1,i:i+n-1);
      end
%
      for k = 1:n,
%
         if (abs(T(i+k-1,i+k-1) - 1) < abs( - T(i+k-1,i+k-1) - 1 ))
           Q(1:m,i+k-1)         = - Q(1:m,i+k-1);
           A(i+k-1,i+k-1:i+n-1) = - A(i+k-1,i+k-1:i+n-1);
           A(i+n:m,i+k-1)  = - A(i+n:m,i+k-1);
           T(i:i+n-1,i+k-1)  = - T(i:i+n-1,i+k-1);
           T(1:i-1,i+k-1) = - T(1:i-1,i+k-1);
           D(i+k-1)      = -1;
         else
           D(i+k-1) = 1;
         end
         T(i+k-1,i+k-1)     = T(i+k-1,i+k-1) - 1;
         T(i+k:i+n-1,i+k-1) = T(i+k:i+n-1,i+k-1) / T(i+k-1,i+k-1);
         A(i+n:m,i+k-1)     = A(i+n:m,i+k-1) / T(i+k-1,i+k-1);
         T(i+k:i+n-1,i+k:i+n-1) = T(i+k:i+n-1,i+k:i+n-1) - T(i+k:i+n-1,i+k-1) * T(i+k-1,i+k:i+n-1);
         A(i+n:m,i+k:i+n-1) = A(i+n:m,i+k:i+n-1) - A(i+n:m,i+k-1) * T(i+k-1,i+k:i+n-1);
      end
%
      for ii=1:n,
         for jj=1:ii-1,
           A(i+ii-1,i+jj-1) = T(i+ii-1,i+jj-1);
         end
      end
%
%    Construction of T
%
      T(1:i-1,i:i+n-1)   = ( - T(1:i-1,i:i+n-1) - T(1:i-1,1:i-1) * A(i:i+n-1,1:i-1)' ) / ( (eye(n,n) + tril(T(i:i+n-1,i:i+n-1),-1) )');
      T(i:i+n-1,i:i+n-1) = - triu( T(i:i+n-1,i:i+n-1) ) / ( (eye(n,n) + tril(A(i:i+n-1,i:i+n-1),-1) )');
%
