function [ A ] = geqrf( A )

   nb = 2;

   [ m, n ] = size( A );

   for j = 1:nb:min( [ m n ]),
   
      jep = j+nb-1;
      if ( jep > n) jep = n; end;
      if ( jep > m) jep = m; end;

      [ A(j:m,j:jep) ] = geqr2( A(j:m,j:jep) );

      [ T ] = larft( A(j:m,j:jep) );

      [ A(j:m,jep+1:n) ] = larfb( A(j:m,j:jep), T, A(j:m,jep+1:n) );

   end

end
