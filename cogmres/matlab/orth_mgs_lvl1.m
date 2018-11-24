%
   function [ Q, r ] = orth_mgs_lvl1( Q )
%
      j = size(Q,2);
      r = zeros(j,1);
%
      for i=1:j-1,
%
         r(i,1) =  Q(:,i)' * Q(:,j);
         Q(:,j) = Q(:,j) - Q(:,i) * r(i,1);
%
      end
%
      r(j,1) = norm( Q(:,j) );
      Q(:,j) = Q(:,j) / r(j,1);
%
   end
%
