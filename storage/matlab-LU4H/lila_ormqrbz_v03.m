%
   function [ B ] = lila_ormqrbz_v03( m, n, k, i, A, ia, ja, lda, B, ib, jb, ldb, T, it, jt, ldt, tb )
%

% tb = [ 4, 5, 3, 5, 7, 9, 5, 4, 5, 7 ];
% nb = [ 20, 13 ];

% tb = [ 4, 5, 3, 5, 7, 9, 5, 4, 5, 7 ];
%
%      [ 3, 1, 5, 3, 1, 4, 3, 4, 9, 5, 4, 5, 7 ];
%        *  ^  ^  ^  ^  *  *  -  -  -  -  -  -

   Q = B;

   ml = m - i - k + 2 ;
 
   ialo = ia+i-1+k-1;
   iahi = m; % or we can write it as ( iahi = ia+i-1+ml-1; ) 
  
   jalo = ja+i-1+k-1;
   jahi = ja+i-1+k-1;
 
   iblo = ib+i-1+k-1;
   ibhi = m; % or we can write it as ( iahi = ia+i-1+ml-1; ) 
  
   jblo = jb;
   jbhi = jb+n-1;
 
   itlo = it+i-1+k-1;
   ithi = it+i-1+k-1;
  
   jtlo = jt+i-1+k-1;
   jthi = jt+i-1+k-1;

  for ii = i+k-1:-1:i,

      V = tril(A(ialo:iahi,jalo:jahi), -1) + eye(size(A(ialo:iahi,jalo:jahi)));

      H = (eye(ml,ml) - V * ( T(itlo:ithi,jtlo:jthi) * V' ) );

      Q(iblo:ibhi,jblo:jbhi) = H*Q(iblo:ibhi,jblo:jbhi);

      iblo = iblo-1;

      ialo = ialo-1;
      jalo = jalo-1;
      jahi = jahi-1;

      itlo = itlo-1;
      ithi = ithi-1;

      jtlo = jtlo-1;
      jthi = jthi-1;

      ml = ml+1;

   end

   B = Q;
%
   end
