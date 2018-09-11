   function [ B ] = lila_ormqrf_v05______( m, n, k, i, A, ia, ja, lda, B, ib, jb, ldb, TT, mt, T, it, jt, ldt, tb )
%
      ialo = ia+i-1;
      iahi = ia+m-1;
      jalo = ja+i-1;
      jahi = jalo+mt-1;
%
      iblo = ib+i-1;
      ibhi = ib+m-1;
      jblo = jb;
      jbhi = jb+n-1;
%
      ilott = 1;
      ihitt = ilott+mt-1;
      jlott = 1;
      jhitt = jlott+mt-1;
%
      ml = m-i+1;
%
      not_done = 1;
%
      jj = mt;
%
       while ( not_done == 1 )

%size(A(ialo:iahi,jalo:jahi))
%ml
%mt
%TT(ilott:ihitt,jlott:jhitt)
fprintf('here\n')
size(eye(ml,ml))
size(TT(ilott:ihitt,jlott:jhitt))
size(A(ialo:iahi,jalo:jahi))
         V = tril(A(ialo:iahi,jalo:jahi), -1) + eye(ml,mt);
size(V)
         H = (eye(ml,ml) - V * ( TT(ilott:ihitt,jlott:jhitt)' * V' ) );
         B(iblo:ibhi,jblo:jbhi) = H*B(iblo:ibhi,jblo:jbhi);

            ialo = ialo + mt;
            iblo = iblo + mt;
            jalo = jalo + mt;
            jlott = jlott + mt;

         if ( jj + mt == k )  % I know the logic will need to change but this is how I manage for now 

            fprintf('here');
            not_done = 0; 

            ml = ml - mt;
            jj = jj + mt;

            jahi = jahi + mt;
            jhitt = jhitt + mt;

            V = tril(A(ialo:iahi,jalo:jahi), -1) + eye(ml,mt);
            H = (eye(ml,ml) - V * ( TT(ilott:ihitt,jlott:jhitt)' * V' ) );
            B(iblo:ibhi,jblo:jbhi) = H*B(iblo:ibhi,jblo:jbhi);

         else if ( jj + mt > k ) % Modular to finish last block for the difference between mt and k
             
             new_mt = mod(k,mt);
             track_t = mt-new_mt;

             ml = ml - new_mt;
             jj = jj + new_mt;

             jahi = jahi+new_mt;
             jhitt = jhitt+new_mt;
             
             ihitt = ihitt-track_t;

         else

             ml = ml - mt;
             jj = jj + mt;

             jahi = jahi + mt;
             jhitt = jhitt + mt;

         end

      end

   end
