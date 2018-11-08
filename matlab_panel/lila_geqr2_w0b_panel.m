%
   function [ A, T, Q, D ] = lila_geqr2_w0b_panel( m, n, i, A, T, Q, D )
%
      work = zeros(i-1,n);
%
      R(1:i-1,i:i+n-1) = A(1:i-1,i:i+n-1);
      for ii=1:n,
         for jj=ii:n,
           R(i+ii-1,i+jj-1) = A(i+ii-1,i+jj-1);
         end
      end
%
      for k = 1:i-1,
         if (D(k) == -1)
           R(k,i:i+n-1) = - R(k,i:i+n-1);
           A(k,i:i+n-1) = - A(k,i:i+n-1);
         end
      end
%
      A(i:m,i:i+n-1)  = Q(i:m,i:i+n-1); 
      work(1:i-1,1:n) = Q(1:i-1,i:i+n-1);
%
      if( i ~= 1),
%
         work(1:i-1,1:n) =  (tril(A(1:i-1,1:i-1),-1) + eye(i-1,i-1)) \ work(1:i-1,1:n);
         A(i:m,i:i+n-1)  = A(i:m,i:i+n-1) - A(i:m,1:i-1) * work(1:i-1,1:n);
%
      end
%
      for k = 1:n,
%
         if (abs(A(i+k-1,i+k-1) - 1) < abs( - A(i+k-1,i+k-1) - 1 ))
%
           Q(1:m,i+k-1)         = - Q(1:m,i+k-1);
           R(i+k-1,i+k-1:i+n-1) = - R(i+k-1,i+k-1:i+n-1);
%
           A(1:m,i+k-1)  = - A(1:m,i+k-1);
           work(1:i-1,k) = - work(1:i-1,k);
           D(i+k-1)      = -1;
%
        else
%
           D(i+k-1) = 1;
%
         end
%
         A(i+k-1,i+k-1)     = A(i+k-1,i+k-1) - 1;
         A(i+k:m,i+k-1)     = A(i+k:m,i+k-1) / A(i+k-1,i+k-1);
         A(i+k:m,i+k:i+n-1) = A(i+k:m,i+k:i+n-1) - A(i+k:m,i+k-1) * A(i+k-1,i+k:i+n-1);
%
      end
%
%    Construction of T
%
      T(1:i-1,i:i+n-1)   = ( - work(1:i-1,1:n) - T(1:i-1,1:i-1) * A(i:i+n-1,1:i-1)' ) / ( (eye(n,n) + tril(A(i:i+n-1,i:i+n-1),-1) )');
      T(i:i+n-1,i:i+n-1) = - triu( A(i:i+n-1,i:i+n-1) ) / ( (eye(n,n) + tril(A(i:i+n-1,i:i+n-1),-1) )');
%
%    Put back the good R
%
      A(1:i-1,i:i+n-1) = R(1:i-1,i:i+n-1);
      for ii=1:n,
         for jj=ii:n,
           A(i+ii-1,i+jj-1) = R(i+ii-1,i+jj-1);
         end
      end
%


% L = tril(A(1:m,1:n),-1)+eye(m,n);
% U = triu(A(1:n,1:n));
% Q - eye(m,n) - ( L * U )
% V = L;
% U = -triu(A(1:n,1:n));
% Q - (  eye(m,n) - V * U )
% T = U / ( V(1:n,1:n) )';
% Q - (  eye(m,n) - V * T * V(1:n,1:n)'  )
% Q - (  eye(m,m) - V * T * V'  )* eye(m,n)
% H= (  eye(m,m) - V * T * V'  ); Q - H* eye(m,n)
