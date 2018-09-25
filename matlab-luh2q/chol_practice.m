%
   clear
%
   n = 11;
%
   A = randn(n);
   A = A*A';
   As = A;
%
%   Method one
%
   for j = 1:n,
      v(j:n,1) = A(j:n,j);
      for k = 1:j-1,
         v(j:n) = v(j:n) - A(j,k)*A(j:n,k);
      end
      A(j:n,j) = v(j:n) / sqrt(v(j));
   end
%
%
%
   fprintf('||GG'' - A|| = %e\n', norm(tril(A)*tril(A)' - As,'fro') / norm(As,'fro'));
%
%   Method two
%
   A = As;
   for j = 1:n,
      if ( j>1 ),
         A(j:n,j) = A(j:n,j) - A(j:n,1:j-1)*A(j,1:j-1)';
      end
      A(j:n,j) = A(j:n,j) / sqrt(A(j,j));
   end
%
%
%
   fprintf('||GG'' - A|| = %e\n', norm(tril(A)*tril(A)' - As,'fro') / norm(As,'fro'));
%
%   Method three (Same as two, without the if statement) 
%
   A = As;
   A(1:n,1) = A(1:n,1) / sqrt(A(1,1));
   for j = 2:n,
      A(j:n,j) = A(j:n,j) - A(j:n,1:j-1)*A(j,1:j-1)';
      A(j:n,j) = A(j:n,j) / sqrt(A(j,j));
   end
%
   fprintf('||GG'' - A|| = %e\n', norm(tril(A)*tril(A)' - As,'fro') / norm(As,'fro'));
