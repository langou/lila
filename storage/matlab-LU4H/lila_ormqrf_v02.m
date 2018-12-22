   function [ B ] = lila_ormqrf_v02( m, n, k, A, ia, ja, lda, B, ib, jb, ldb, T, it, jt, ldt, ub )
%
      iblo = ib;
      ibhi = ib+m-1;
      jblo = jb;
      jbhi = jb+n-1;
%
      ialo = ia;
      iahi = ia+m-1;
      jalo = ja;
%
      itlo = it;
      jtlo = jt;
%
      ml =m;
%
%     ib = [ 4, 7, 9, 3, 5, 20 ] this is the ib fed as ub
      ub = [ 4, 7, 5, 4, 3, 3 ]; % Breaking 9 into 5,4

      

      p = 1;
      while( jtlo + ub(p) - 1 <= jt )
         p = p + 1;
         jahi = jalo+ub(p)-1
         ithi = itlo+ub(p)-1
         jthi = jtlo+ub(p)-1
         ialo = ialo+ub(p)
         jalo = jalo+ub(p)
         iblo = iblo+ub(p)
         itlo = itlo+ub(p)
         jtlo = jtlo+ub(p)
         ml = ml - ub(p)
      end

%
      for i = 1:3,
%
      jahi = jalo+ub(i)-1;
      ithi = itlo+ub(i)-1;
      jthi = jtlo+ub(i)-1;
%
         V = tril(A(ialo:iahi,jalo:jahi), -1) + eye(ml,ub(i));
         H = (eye(ml,ml) - V * ( T(itlo:ithi,jtlo:jthi)' * V' ) );
         B(iblo:ibhi,jblo:jbhi) = H*B(iblo:ibhi,jblo:jbhi);
%
      ialo = ialo+ub(i);
      jalo = jalo+ub(i);
      iblo = iblo+ub(i);
      itlo = itlo+ub(i);
      jtlo = jtlo+ub(i);
      ml = ml - ub(i);
%
      end
%











   end
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      jahi = jalo+ub(1)-1;
%      ithi = itlo+ub(1)-1;
%      jthi = jtlo+ub(1)-1;
%      ialo = ialo+ub(1);
%      jalo = jalo+ub(1);
%      iblo = iblo+ub(1);
%      itlo = itlo+ub(1);
%      jtlo = jtlo+ub(1);
%      ml = ml - ub(1);

