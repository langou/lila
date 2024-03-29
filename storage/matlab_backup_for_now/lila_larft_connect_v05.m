   function [ TT ] = lila_larft_connect_v05( m, n, i, A, TT, mt )
%
   vb= mt - mod(i-1,mt) ;
   if ( vb > n ), vb = n; end
%
   itlo = mod(i-1,mt)+1;
   ithi = itlo+vb-1;
   jtlo = i;
   jthi = jtlo+vb-1;
%
   TT(1:itlo-1,jtlo:jthi) = ...
       - ( TT(1:itlo-1,jtlo-(itlo)+1:jtlo-1) * A(jtlo:m,jtlo-(itlo)+1:jtlo-1)' )* ...
           ( (tril(A(jtlo:m,jtlo:jthi),-1) + eye(size((A(jtlo:m,jtlo:jthi)))) ) * TT(itlo:ithi,jtlo:jthi) );
%
   end
