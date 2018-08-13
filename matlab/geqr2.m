function [ A ] = geqr2( A )
   [ m, n ] = size( A );
   for j = 1:min( [ m n ]),
   
      [ A(j:m,j) ] = larfg( A(j:m,j) );

      if (j < n ), [ A(j:m,j+1:n) ] = larfL( A(j:m,j), A(j:m,j+1:n) ); end

   end
end
