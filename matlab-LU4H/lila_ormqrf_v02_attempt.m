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
%%%%
%
%  working through some ideas
%
%
      ll = size(ub,2);
      count = 0;
      kill_num = 0;
      tmp_val1 = 0;
%
      for z = 1:ll,
         if tmp_val1 < k, tmp_val1 = tmp_val1 + ub(z); count = count+1; else end;
         if tmp_val1 > k, tmp_val1 = tmp_val1 - ub(z); count = count-1; else end;
      end
%
      for z = 1:count,
         kill_num = kill_num + ub(z);
      end
      tmp_val2 = abs(k - kill_num); % tmp_val2 is the value to move backwards in ib for the correct block size
%
      jahi = ja+ub(1)-1;
      jthi = jt+ub(1)-1;
      ithi = it+ub(1)-1;
      jbhi = jb+ub(1)-1;

      for z = 1:count-1,
      
%         V = tril(A(ialo:iahi,jalo:jahi), -1) + eye(m,k);
         V = tril(A(ialo:iahi,jalo:jahi), -1) + eye(m,ub(z));
%
T(itlo:ithi,jtlo:jthi)'
%
         H = (eye(m,m) - V(ialo:iahi,jalo:jahi) * ( T(itlo:ithi,jtlo:jthi)' * V(ialo:iahi,jalo:jahi)' ) );

%         B(iblo:ibhi,jblo:jbhi) = H*B(iblo:ibhi,jblo:jbhi);
         B(iblo:ibhi,jblo:jbhi) = H*B(iblo:ibhi,jblo:jbhi);

         jahi = ub(z+1);
         jtlo = jtlo + ub(z);
         itlo = itlo + ub(z);
         jthi = jthi + ub(z+1);
         ithi = ithi + ub(z+1);
         jblo = jb + ub(z);
         jbhi = jb+ub(z+1);

      end
%      
      if tmp_val2 == 0,
         % do nothing - we've exhausted the block
      else
%
         jahi = jahi - tmp_val2;
         jthi = jthi - tmp_val2;
         ithi = ithi - tmp_val2;
         jbhi = jb - tmp_val2;

%
         V = tril(A(ialo:iahi,jalo:jahi), -1) + eye(m,ub(count) - tmp_val2);
%
         H = (eye(m,m) - V(ialo:iahi,jalo:jahi) * ( T(itlo:ithi,jtlo:jthi)' * V(ialo:iahi,jalo:jahi)' ) );
%
%         B(iblo:ibhi,jblo:jbhi) = H*B(iblo:ibhi,jblo:jbhi);
         B(iblo:ibhi,jblo:jbhi) = H*B(iblo:ibhi,jblo:jbhi);
%
       end
%
T(itlo:ithi,jtlo:jthi)'
%
   end

