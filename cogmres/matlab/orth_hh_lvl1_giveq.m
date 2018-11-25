%
   function [ a ] = orth_hh_lvl1_giveq( V, tau )
%
      m = size(V,1);
      j = size(V,2);
%
      a = zeros(m,1);
%
      a(j,1) = - tau(j) ;
      a(j+1:m,1) = V(j+1:m,j) * a(j,1) ;
      a(j,1) = 1.0e+00 + a(j,1) ;
      for i = j-1:-1:1,
        a(i,1) = - tau(i) * ( V(i+1:m,i)' * a(i+1:m,1) ) ;
        a(i+1:m,1) = a(i+1:m,1) + V(i+1:m,i) * a(i,1) ;
      end
%
   end
%

