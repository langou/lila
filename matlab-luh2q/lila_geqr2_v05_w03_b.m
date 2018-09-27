%
   function [ A, T, Q ] = lila_geqr2_v05_w03_b( m, n, i, mt, A, T, Q )
%
   Q(i:m,i:i+n-1) = A(i:m,i:i+n-1);
%
%     Cholesky QR
%
   R(1:n,1:n) = Q(i:m,i:i+n-1)'*Q(i:m,i:i+n-1); 
   R(1:n,1:n) = chol( R(1:n,1:n), 'upper' ); 
   Q(i:m,i:i+n-1) = Q(i:m,i:i+n-1) / R(1:n,1:n);
%
   R2(1:n,1:n) = Q(i:m,i:i+n-1)'*Q(i:m,i:i+n-1); 
   R2(1:n,1:n) = chol( R2(1:n,1:n), 'upper' ); 
   Q(i:m,i:i+n-1) = Q(i:m,i:i+n-1) / R2(1:n,1:n);
   R(1:n,1:n) = R2(1:n,1:n) * R(1:n,1:n);
   clear R2
%
%     Householder reconstruction
%
   A(i:m,i:i+n-1) = Q(i:m,i:i+n-1); 
   D = -eye(n,n);
%
   for k = 1:n,
      if (abs(1 - A(k,k)) < abs( 1 + A(k,k) ))
         Q(i:m,i+k-1) = -Q(i:m,i+k-1);
         R(k,k:n) = -R(k,k:n);
      else 
         A(i:m,i+k-1) = -A(i:m,i+k-1);
         D(k,k) = 1.0e+00;
      end
      A(i+k-1,i+k-1) = 1 + A(i+k-1,i+k-1);
      A(i+k:m,i+k-1) = A(i+k:m,i+k-1) / A(i+k-1,i+k-1);
      A(i+k:m,i+k:i+n-1) = A(i+k:m,i+k:i+n-1) - A(i+k:m,i+k-1) * A(i+k-1,i+k:i+n-1);
   end
%
%    Construction of T with the structure of mt
%
   vb = mt - mod(i-1,mt) ;
%
   if ( vb > n ), vb = n; end
%
   itlo = mod(i-1,mt)+1; if (itlo == 0), itlo = mt, end;
   ithi = itlo+vb-1;
   jtlo = i;
   jthi = jtlo+vb-1;
%
   T(itlo:ithi,jtlo:jthi) = triu( A(jtlo:jthi,jtlo:jthi) ) / ( (eye(vb,vb) + tril(A(jtlo:jthi,jtlo:jthi),-1) )');
%
   itlo = mod(itlo+vb-1,mt)+1;
   jtlo = jtlo+vb;
%
   if( jtlo+mt-1 >= i+n-1 ) vb = i+n-1+1-jtlo; else, vb = mt; end
   ithi = itlo+vb-1;
   jthi = jtlo+vb-1;
%
   while (( itlo == 1 )&&(vb~=0))
%
      T(itlo:ithi,jtlo:jthi) = triu( A(jtlo:jthi,jtlo:jthi) ) / ( (eye(vb,vb) + tril(A(jtlo:jthi,jtlo:jthi),-1) )');
%
      itlo = mod(itlo+vb-1,mt)+1; if (itlo == 0), itlo = mt, end;
      jtlo = jtlo+vb;
      if( jtlo+mt-1 >= i+n-1) vb = i+n-1+1-jtlo; else, vb = mt; end
      ithi = itlo+vb-1;
      jthi = jtlo+vb-1;
%
   end
%
%    Put back the good R
%
   for ii=1:n,
      for jj=ii:n,
        A(i+ii-1,i+jj-1) = R(ii,jj);
      end
   end
%
   end
