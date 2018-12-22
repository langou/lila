%
   function [ B ] = lapack_larfL( a, B )
%
      [m,n]=size(B);
%
      tau = 2 / ( 1 + norm(a(2:m,1))^2 );
%
      tmp(1,1:n) = B(1,1:n) + a(2:m,1)' * B(2:m,1:n);
      tmp(1,1:n) = tau * tmp(1,1:n);
      B(1,1:n) = B(1,1:n) - tmp(1,1:n);
      B(2:m,1:n) = B(2:m,1:n) - a(2:m,1) * tmp(1,1:n);
%
   end
