%
   function [ q, r ] = orth_cgs( Q, q )
%
      j = size(Q,2)+1;
      r = zeros(j,1);
%
      r(1:j-1,1) = Q(:,1:j-1)' * q;
      q = q - Q(:,1:j-1) * r(1:j-1,1);
%
      r(j,1) = norm( q );
      q = q / r(j,1);
%
   end
%
