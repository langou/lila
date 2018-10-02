%
   function [ Q ] = lila_ormqrbz_v05_w01( m, n, k, i, j, mt, A, T, Q )
%
%      Q(i:i+k-1,j:j+n-1) = zeros( size( Q(i:i+k-1,j:j+n-1) ));
%      for ii = i+k-1:-1:i,
%         V = tril(A(ii:m,ii), -1) + eye(size(A(ii:m,1)));
%         H = (eye(m-ii+1,m-ii+1) - V * ( T(1,ii) * V' ) );
%         Q(ii:m,j:j+n-1) = H * Q(ii:m,j:j+n-1);
%      end
%
%      Q(i:i+k-1,j:j+n-1) = zeros( size( Q(i:i+k-1,j:j+n-1) ));
%      for ii = i+k-1:-1:i,
%         V = tril(A(ii:m,ii), -1) + eye(size(A(ii:m,1)));
%         Q(ii:m,j:j+n-1) = (eye(m-ii+1,m-ii+1) - V * ( T(1,ii) * V' ) ) * Q(ii:m,j:j+n-1);
%      end
%
%     Q(i:i+k-1,j:j+n-1) = zeros( size( Q(i:i+k-1,j:j+n-1) ));
%     for ii = i+k-1:-1:i,
%        Q(ii:m,j:j+n-1) = Q(ii:m,j:j+n-1) - (tril(A(ii:m,ii), -1) + eye(size(A(ii:m,1)))) * ( T(1,ii) * ( (tril(A(ii:m,ii), -1) + eye(size(A(ii:m,1))))' * Q(ii:m,j:j+n-1) ) );
%     end
%
%     work(1,1:n) = zeros(1,n);
%     Q(i:i+k-1,j:j+n-1) = zeros( size( Q(i:i+k-1,j:j+n-1) ));
%     for ii = i+k-1:-1:i,
%        work(1,1:n) = [ 1; A(ii+1:m,ii) ]' * Q(ii:m,j:j+n-1);
%        work(1,1:n) = T(1,ii) * work(1,1:n);
%        Q(ii:m,j:j+n-1) = Q(ii:m,j:j+n-1) - [ 1; A(ii+1:m,ii) ] * work(1,1:n);
%     end
%
%     work(1,1:n) = zeros(1,n);
%     for ii = i+k-1:-1:i,
%        work(1,1:n) = [ A(ii+1:m,ii) ]' * Q(ii+1:m,j:j+n-1);
%        work(1,1:n) = T(1,ii) * work(1,1:n);
%        Q(ii,j:j+n-1) = - work(1,1:n);
%        Q(ii+1:m,j:j+n-1) = Q(ii+1:m,j:j+n-1) - [ A(ii+1:m,ii) ] * work(1,1:n);
%     end
%
      for ii = i+k-1:-1:i,
         Q(ii,j:j+n-1) = A(ii+1:m,ii)' * Q(ii+1:m,j:j+n-1);
         Q(ii,j:j+n-1) = - T(1,ii) * Q(ii,j:j+n-1);
         Q(ii+1:m,j:j+n-1) = Q(ii+1:m,j:j+n-1) + A(ii+1:m,ii) * Q(ii,j:j+n-1);
      end
%
   end
