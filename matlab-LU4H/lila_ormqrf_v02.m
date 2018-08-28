   function [ B ] = lila_ormqrf_v02( m, n, k, A, ia, ja, lda, B, ib, jb, ldb, T, it, jt, ldt, ub, nb, l )
%
%  Because you have ib already I've changed it to ub for user block
%  Added nb to loop through correctly
%  and the loop variable k - but was taken so named it l
%
      ialo = ia;
      iahi = ia+m-1;
      jalo = ja;
      jahi = ja+k-1;
%
      iblo = ib;
      ibhi = ib+m-1;
      jblo = jb;
      jbhi = jb+n-1;
%
      itlo = it;
      ithi = it+k-1;
      jtlo = jt;
      jthi = jt+k-1;
%
      

%      V = tril(A(ialo:iahi,jalo:jahi), -1) + eye(m,k);
      V = tril(A(ialo:iahi,), -1) + eye(m,k);

      H = (eye(m,m) - V * ( T(itlo:ithi,jtlo:jthi)' * V' ) );
%      H = (eye(m,m) - V * ( T()' * V' ) );

      B(iblo:ibhi,jblo:jbhi) = H*B(iblo:ibhi,jblo:jbhi);
%
   end
