%
   function [ A, T, Q, D ] = lila_geqr2_w0b_panel_mt( m, n, i, j, mt, A, T, Q, D )
%
%     lila_ormhr_w0b % ormqr
      work = zeros( i-1, n);
      zork = zeros( j-i+1-1, n);
%
      itlo = mod(i-1,mt)+1; if (itlo == 0), itlo = mt, end;
%
      for k = 1:j-1,
         if (D(k) == -1)
           A(k,j:j+n-1) = - A(k,j:j+n-1);
         end
      end
%
%      work(1:i-1,1:n)    = Q(1:i-1,j:j+n-1);
      zork(1:j-i+1-1,1:n)= Q(i:j-1,j:j+n-1);
      T(itlo:itlo+n-1,j:j+n-1) = Q(j:j+n-1,j:j+n-1);
      A(j+n:m,j:j+n-1)   = Q(j+n:m,j:j+n-1); 
%
%      work(1:i-1,1:n)    = (tril(A(1:i-1,1:i-1),-1) + eye(i-1,i-1)) \ work(1:i-1,1:n);
%      zork(1:j-i+1-1,1:n)   = zork(1:j-i+1-1,1:n) - A(i:j-1,1:i-1) * work(1:i-1,1:n);
      zork(1:j-i+1-1,1:n)   = (tril(A(i:j-1,i:j-1),-1) + eye(j-i,j-i)) \ zork(1:j-i+1-1,1:n);
%      T(j:j+n-1,j:j+n-1) = T(j:j+n-1,j:j+n-1) - A(j:j+n-1,1:i-1) * work(1:i-1,1:n);
      T(itlo:itlo+n-1,j:j+n-1) = T(itlo:itlo+n-1,j:j+n-1) - A(j:j+n-1,i:j-1) * zork(1:j-i+1-1,1:n);
%      A(j+n:m,j:j+n-1)   = A(j+n:m,j:j+n-1) - A(j+n:m,1:i-1) * work(1:i-1,1:n);
      A(j+n:m,j:j+n-1)   = A(j+n:m,j:j+n-1) - A(j+n:m,i:j-1) * zork(1:j-i+1-1,1:n);
%
%
%     lila_orgh2_w0b_panel
      for k = 1:n,
%
      if (abs(T(itlo+k-1,j+k-1) - 1) < abs( - T(itlo+k-1,j+k-1) - 1 ))
         T(itlo:itlo+n-1,j+k-1)     = - T(itlo:itlo+n-1,j+k-1);
         zork(1:j-i+1-1,k)   = - zork(1:j-i+1-1,k);
         D(j+k-1)             = -1;
      else
         D(j+k-1)             = 1;
      end
         T(itlo+k-1,j+k-1)         = T(itlo+k-1,j+k-1) - 1;
         T(itlo+k:itlo+n-1,j+k-1)     = T(itlo+k:itlo+n-1,j+k-1) / T(itlo+k-1,j+k-1);
         T(itlo+k:itlo+n-1,j+k:j+n-1) = T(itlo+k:itlo+n-1,j+k:j+n-1) - T(itlo+k:itlo+n-1,j+k-1) * T(itlo+k-1,j+k:j+n-1);
      end
%
      for k = 1:n,
         if ( D(j+k-1) == -1),
            A(j+k-1,j+k-1:j+n-1) = - A(j+k-1,j+k-1:j+n-1);
         end
      end
%
      for k = 1:n,
         if ( D(j+k-1) == -1),
            A(j+n:m,j+k-1) = - A(j+n:m,j+k-1);
         end
      end
%
      for k = 1:n,
         if ( D(j+k-1) == -1),
            Q(i:m,j+k-1) = - Q(i:m,j+k-1);
         end
      end
%
%     Putting back the lower triangle of V back into A
%
      for ii=1:n,
         for jj=1:ii-1,
           A(j+ii-1,j+jj-1) = T(itlo+ii-1,j+jj-1);
         end
      end

     A(j+n:m,j:j+n-1) = A(j+n:m,j:j+n-1) / triu( T(itlo:itlo+n-1,j:j+n-1) );

%
      T(itlo:itlo+n-1,j:j+n-1) = - triu( T(itlo:itlo+n-1,j:j+n-1) ) / ( (eye(n,n) + tril(A(j:j+n-1,j:j+n-1),-1) )');
%
