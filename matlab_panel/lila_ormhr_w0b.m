%
   function [ A, T ] = lila_ormhr_w0b( m, n, i, j, mt, A, T, Q, D )
%
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
   end
