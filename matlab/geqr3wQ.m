function [ Q, R, V, T ] = geqr3wQ ( A )

   [ m, n ] = size(A);

   [ A ] = geqrf( A );

   V = tril( A, -1 ) + eye( m, n );
   T = larft( V );
   R = triu(A(1:n,1:n));

   Q = eye(m,n);
   [ Q ] = larfb( V, T', Q );

end
