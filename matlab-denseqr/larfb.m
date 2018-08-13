function [ A ] = larfb( V, T, A )

   [ m, n] = size(A);
   [ m, k] = size(V);

   for j = 1:k,
      for i = 1:j-1,
         V(i,j) = 0.0e+00;
      end
      V(j,j) = 1.0e+00;
   end

   A = A - V * ( T' * ( V' * A ) ) ;

end

