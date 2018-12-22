function [ B ] = larfR( a, B )

   [m,n]=size(B);

   tau = 2 / ( 1 + norm(a(2:n,1))^2 );
   tmp = B(1:m,1) + B(1:m,2:n) * a(2:n,1) ;
   tmp = tau * tmp ;
   B(1:m,1) = B(1:m,1)  - tmp ;
   B(1:m,2:n) = B(1:m,2:n)  - tmp * a(2:n,1)' ;

end
