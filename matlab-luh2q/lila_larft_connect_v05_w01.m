   function [ T ] = lila_larft_connect_v05_w01( m, n, i, mt, A, T )
%
%
   itlo = i;
   ithi = i;
   jtlo = i;
   jthi = j;
%
   for ii = 1:n,
%
   T(1:itlo-1,jtlo:jthi) = ...
       - ( T(1:itlo-1,jtlo-(itlo)+1:jtlo-1) * A(jtlo:m,jtlo-(itlo)+1:jtlo-1)' )* ...
           ( (tril(A(jtlo:m,jtlo:jthi),-1) + eye(size((A(jtlo:m,jtlo:jthi)))) ) * T(itlo:ithi,jtlo:jthi) );
%
   end
%
   end
