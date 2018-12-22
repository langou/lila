function [ T ] = lapack_larft( V )

   [ m, k] = size(V);

   for j = 1:k,
      for i = 1:j-1,
         V(i,j) = 0.0e+00;
      end
      V(j,j) = 1.0e+00;
   end

   T = zeros( k, k);

   for j = 1:k,

      tau(j) = 2 / ( 1 + norm(V(j+1:m,j))^2 );
      T(1:j-1,j) = -T(1:j-1,1:j-1)*V(1:m,1:j-1)'*V(1:m,j)*tau(j);
      T(j,j) = tau(j);

   end


end
