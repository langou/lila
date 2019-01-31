%
   function [ A, T, Q, D ] = lila_geqr2_w0b_panel_mt( m, n, i, j, mt, A, T, Q, D )
%
%     lila_ormhr_w0b % ormqr
      work = zeros( i-1, n);
%
      for k = 1:j-1,
         if (D(k) == -1)
           A(k,j:j+n-1) = - A(k,j:j+n-1);
         end
      end
%
      work(1:i-1,1:n)    = Q(1:i-1,j:j+n-1);
      T(i:j-1,j:j+n-1)   = Q(i:j-1,j:j+n-1);
      T(j:j+n-1,j:j+n-1) = Q(j:j+n-1,j:j+n-1);
      A(j+n:m,j:j+n-1)   = Q(j+n:m,j:j+n-1); 
%
      work(1:i-1,1:n)    = (tril(A(1:i-1,1:i-1),-1) + eye(i-1,i-1)) \ work(1:i-1,1:n);
      T(i:j-1,j:j+n-1)   = T(i:j-1,j:j+n-1) - A(i:j-1,1:i-1) * work(1:i-1,1:n);
      T(i:j-1,j:j+n-1)   = (tril(A(i:j-1,i:j-1),-1) + eye(j-i,j-i)) \ T(i:j-1,j:j+n-1);
      T(j:j+n-1,j:j+n-1) = T(j:j+n-1,j:j+n-1) - A(j:j+n-1,1:i-1) * work(1:i-1,1:n);
      T(j:j+n-1,j:j+n-1) = T(j:j+n-1,j:j+n-1) - A(j:j+n-1,i:j-1) * T(i:j-1,j:j+n-1);
      A(j+n:m,j:j+n-1)   = A(j+n:m,j:j+n-1) - A(j+n:m,1:i-1) * work(1:i-1,1:n);
      A(j+n:m,j:j+n-1)   = A(j+n:m,j:j+n-1) - A(j+n:m,i:j-1) * T(i:j-1,j:j+n-1);
%
%
%     lila_orgh2_w0b_panel
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
%
%     lila_larft_w0b_connect
      T(i:j-1,j:j+n-1)   = ( - T(i:j-1,j:j+n-1) - T(i:j-1,1:j-1) * A(j:j+n-1,1:j-1)' ) ;
      T(i:j-1,j:j+n-1)   = T(i:j-1,j:j+n-1) / ( (eye(n,n) + tril(A(j:j+n-1,j:j+n-1),-1) )');
