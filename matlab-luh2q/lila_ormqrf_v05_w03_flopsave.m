%
   function [ A ] = lila_ormqrf_v05_w03_flopsave( m, n, k, i, j, mt, A, T )
%
      vb= mt - mod(i-1,mt) ;
      if ( vb > k ), vb = k; end
%
      ialo = i;
      iahi = m;
      jalo = i;
      jahi = jalo+vb-1;
%
      iblo = i;
      ibhi = m;
      jblo = j;
      jbhi = j+n-1;
%
      itlo = mod(i-1,mt)+1;
      ithi = itlo+vb-1;
      jtlo = i;
      jthi = jtlo+vb-1;
%
      ml = m-i+1;
%
      not_done = 1;
%
      jj = 1;
%
      while ( not_done == 1 )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%         V = tril(A(ialo:iahi,jalo:jahi), -1) + eye(ml,vb);
%         H = (eye(ml,ml) - V * ( T(itlo:ithi,jtlo:jthi)' * V' ) );
%         A(iblo:ibhi,jblo:jbhi) = H*A(iblo:ibhi,jblo:jbhi);

%         V = tril(A(ialo:iahi,jalo:jahi), -1) + eye(ml,vb);
%         A(iblo:ibhi,jblo:jbhi) = ( A(iblo:ibhi,jblo:jbhi) -  V * ( T(itlo:ithi,jtlo:jthi)' * V' ) * A(iblo:ibhi,jblo:jbhi));

%         V = tril(A(ialo:iahi,jalo:jahi), -1) + eye(ml,vb);
%         A(iblo:ibhi,jblo:jbhi) = ( A(iblo:ibhi,jblo:jbhi) -  (tril(A(ialo:iahi,jalo:jahi), -1) + eye(ml,vb)) * ( T(itlo:ithi,jtlo:jthi)' * (tril(A(ialo:iahi,jalo:jahi), -1) + eye(ml,vb))' ) * A(iblo:ibhi,jblo:jbhi));
       
%         V = tril(A(ialo:iahi,jalo:jahi), -1) + eye(ml,vb);
%         work =  V' * A(iblo:ibhi,jblo:jbhi);
%         work =  T(itlo:ithi,jtlo:jthi)' * work;
%         A(iblo:ibhi,jblo:jbhi) = A(iblo:ibhi,jblo:jbhi) - V * work;

%         work =  (tril(A(ialo:iahi,jalo:jahi), -1) + eye(ml,vb))' * A(iblo:ibhi,jblo:jbhi);
%         work =  T(itlo:ithi,jtlo:jthi)' * work;
%         A(iblo:ibhi,jblo:jbhi) = A(iblo:ibhi,jblo:jbhi) - (tril(A(ialo:iahi,jalo:jahi), -1) + eye(ml,vb)) * work;

%         work(1:vb,1:n) =  (tril(A(jalo:jahi,jalo:jahi), -1) + eye(vb,vb))' * A(iblo:iblo+vb-1,jblo:jbhi);
%         work(1:vb,1:n) =  work(1:vb,1:n) + A(jahi+1:iahi,jalo:jahi)' * A(iblo+vb:ibhi,jblo:jbhi);
%         work(1:vb,1:n) =  T(itlo:ithi,jtlo:jthi)' * work(1:vb,1:n);
%         A(iblo:ibhi,jblo:jbhi) = A(iblo:ibhi,jblo:jbhi) - (tril(A(ialo:iahi,jalo:jahi), -1) + eye(ml,vb)) * work(1:vb,1:n);

%         work(1:vb,1:n)              =  (tril(A(jalo:jahi,jalo:jahi), -1) + eye(vb,vb))' * A(iblo:iblo+vb-1,jblo:jbhi);
%         work(1:vb,1:n)              =  work(1:vb,1:n) + A(jahi+1:iahi,jalo:jahi)' * A(iblo+vb:ibhi,jblo:jbhi);
%         work(1:vb,1:n)              =  T(itlo:ithi,jtlo:jthi)' * work(1:vb,1:n);
%         A(iblo:iblo+vb-1,jblo:jbhi) =  A(iblo:iblo+vb-1,jblo:jbhi) - (tril(A(ialo:jahi,jalo:jahi), -1) + eye(vb,vb))*work(1:vb,1:n);
%         A(iblo+vb:ibhi,jblo:jbhi)   =  A(iblo+vb:ibhi,jblo:jbhi) - A(jahi+1:iahi,jalo:jahi) * work(1:vb,1:n);
%
         lda = -1;
         ldt = -1;
         ldwork = -1;
%
        work = A(iblo:iblo+vb-1,jblo:jbhi);
        [ work ] = blas_trmm( 'L', 'L', 'T', 'U', vb, n, (+1.0e00), A, jalo, jalo, lda, work, 1, 1, ldwork );
        [ work ] = blas_gemm( 'T', 'N', vb, n, m-jahi, (+1.0e00), A, jahi+1, jalo, lda, A, iblo+vb, jblo, lda, (+1.00e00), work, 1, 1, ldwork );
        [ work ] = blas_trmm( 'L', 'U', 'T', 'N', vb, n, (+1.0e00), T, itlo, jtlo, ldt, work, 1, 1, ldwork );
        [ A ] = blas_gemm( 'N', 'N', m-vb-iblo+1, n, vb, (-1.0e00), A, jahi+1, jalo, lda, work, 1, 1, ldwork, (+1.00e00), A, iblo+vb, jblo, lda );
        [ work ] = blas_trmm( 'L', 'L', 'N', 'U', vb, n, (-1.0e00), A, ialo, jalo, lda, work, 1, 1, ldwork );
        A(iblo:iblo+vb-1,jblo:jbhi) =  A(iblo:iblo+vb-1,jblo:jbhi) + work;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

         if ( jj + vb - 1 == k ) 

            not_done = 0; 

         else if ( jj + vb - 1 > k ) 

             fprintf('defensive programming, we should never have been there, abort\n'); return;

         else

            ml = ml-vb;

            ialo = ialo+vb;
            iblo = iblo+vb;
            jalo = jalo+vb;
            itlo = 1;
            jtlo = jtlo+vb;

            jj = jj + vb;
 
            if ( ( jj+mt-1 ) <= k ), vb = mt; else vb = k-jj+1; end;

            jahi = jahi+vb;
            ithi = vb;
            jthi = jthi+vb;

         end
%
      end
%
   end
