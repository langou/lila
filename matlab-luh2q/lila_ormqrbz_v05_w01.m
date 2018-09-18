%
   function [ Q ] = lila_ormqrbz_v05_w01( m, n, k, i, jb, mt, A, T, Q )
%
      iblo = i;
%
      jblo = jb;
      jbhi = jb+n-1;
%
      Q(iblo:i+k-1,jblo:jbhi) = zeros( size( Q(i:i+k-1,jblo:jbhi) ));

%
      ml = m-i-k+2;
% 
      ialo = i+k-1;
      iahi = m;  
%  
      jalo = i+k-1;
      jahi = i+k-1;
% 
      iblo = i+k-1;
      ibhi = m; 
%  
      jblo = jb;
      jbhi = jb+n-1;
% 
      itlo = i+k-1;
      ithi = i+k-1;
%  
      jtlo = i+k-1;
      jthi = i+k-1;
%
      for ii = k:-1:1,
%
         V = tril(A(ialo:iahi,jalo:jahi), -1) + eye(size(A(ialo:iahi,jalo:jahi)));
%
         H = (eye(ml,ml) - V * ( lapack_larft(V) * V' ) );
%
         Q(iblo:ibhi,jblo:jbhi) = H * Q(iblo:ibhi,jblo:jbhi);
%
         iblo = iblo-1;
%
         ialo = ialo-1;
         jalo = jalo-1;
         jahi = jahi-1;
%
         itlo = itlo-1;
         ithi = ithi-1;
%
         jtlo = jtlo-1;
         jthi = jthi-1;
%
         ml = ml+1;
% 
      end
%
   end
