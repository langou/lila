%
   function [ A, T, Q ] = lila_geqr2_v05_q00( m, n, i, mt, A, T, Q )
%
   As = A;
%  [ A ] = lapack_geqr2( m, n, i, A ); 
%
%  [ T ] = lila_larft_v05_w03( m, n, i, mt, A, T );  
%
%  [ Q ] = lila_orgqrf_v05_w03( m, n, i, mt, A, T, Q );  
%
%
   ilo = i;
   ihi = i+n-1;
   ml = m-i+1;
%
%  [ Q(i:m,ilo:ihi), R ] = qr( As(i:m,ilo:ihi), 0 );

   R(1:n,1:n) = As(i:m,ilo:ihi)'*As(i:m,ilo:ihi); 
   R(1:n,1:n) = chol( R(1:n,1:n), 'upper' ); 
   Q(i:m,ilo:ihi) = As(i:m,ilo:ihi) / R(1:n,1:n);
%
   R2(1:n,1:n) = Q(i:m,ilo:ihi)'*Q(i:m,ilo:ihi); 
   R2(1:n,1:n) = chol( R2(1:n,1:n), 'upper' ); 
   Q(i:m,ilo:ihi) = Q(i:m,ilo:ihi) / R2(1:n,1:n);
   R(1:n,1:n) = R2(1:n,1:n) * R(1:n,1:n);
%



   QQQ = Q(i:m,ilo:ihi); Rs = R(1:n,1:n);
   QQQs = QQQ;
%
   mmm = m-i+1;
   nnn = n;
   D = -eye(nnn,nnn);
   for k = 1:nnn,
      if (abs(1 - QQQ(k,k)) < abs( 1 + QQQ(k,k) )) 
         Q(i:m,i+k-1) = -Q(i:m,i+k-1);
         R(k,k:nnn) = -R(k,k:nnn);
      else 

         QQQ(1:mmm,k) = -QQQ(1:mmm,k);
         D(k,k) = 1.0e+00;
         %Q(i:m,i+k-1) = -Q(i:m,i+k-1);
         %R(k,k:nnn) = -R(k,k:nnn);

      end
      QQQ(k,k) = 1 + QQQ(k,k);
      QQQ(k+1:mmm,k) = QQQ(k+1:mmm,k) / QQQ(k,k);
      QQQ(k+1:mmm,k+1:nnn) = QQQ(k+1:mmm,k+1:nnn) -  QQQ(k+1:mmm,k) * QQQ(k,k+1:nnn);
   end
%
%  Checks and construction of H
%
   U = triu(QQQ); U=U(1:nnn,1:nnn);
   L = tril(QQQ,-1)+eye(mmm,nnn);
%  fprintf('|| QQQ - L*U ||             = %d \n', norm( ( eye(mmm,nnn) - QQQs * D ) - L*U,'fro') );
%
   V = L;
   TTT = U / (V(1:nnn,1:nnn)');


   for ii=1:n,
      for jj=ii:n,
        A(i+ii-1,i+jj-1) = R(ii,jj);
      end
   end



   for ii=1:nnn,
      for jj=1:ii-1,
        A(i+ii-1,i+jj-1) = QQQ(ii,jj);
      end
   end

   for ii=nnn+1:mmm,
      for jj=1:nnn,
        A(i+ii-1,i+jj-1) = QQQ(ii,jj);
      end
   end


   [ T ] = lila_larft_v05_w03( m, n, i, mt, A, T );  

%
   end
