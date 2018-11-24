%
   function [ q, r ] = orth_mgs_lvl1( Q, q )
%
      j = size(Q,2)+1;
      r = zeros(j,1);
%
      for i=1:j-1,
%
         r(i,1) =  Q(:,i)' * q;
         q = q - Q(:,i) * r(i,1);
%
      end
%
      r(j,1) = norm( q );
      q = q / r(j,1);
%
   end
%
