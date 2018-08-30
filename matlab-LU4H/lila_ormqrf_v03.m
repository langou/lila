   function [ B ] = lila_ormqrf_v03( m, n, k, i, A, ia, ja, lda, B, ib, jb, ldb, T, it, jt, ldt, tb )

%%%%%%%%%%%%%%%%%%%

%     ialo = ia+i-1;
%     iahi = ia+m-1;
%     jalo = ja+i-1;
%     jahi = jalo+k-1;
%
%     iblo = ib+i-1;
%     ibhi = ib+m-1;
%     jblo = jb;
%     jbhi = jb+n-1;
%
%     itlo = it+i-1;
%     ithi = itlo+k-1;
%     jtlo = jt+i-1;
%     jthi = jtlo+k-1;
%
%     ml = m-i+1;
%
%     V = tril(A(ialo:iahi,jalo:jahi), -1) + eye(ml,k);
%     H = (eye(ml,ml) - V * ( T(itlo:ithi,jtlo:jthi)' * V' ) );
%     B(iblo:ibhi,jblo:jbhi) = H*B(iblo:ibhi,jblo:jbhi);

%%%%%%%%%%%%%%%%%%%

%     ib1 = 1;
%     ialo = ia+i-1;
%     iahi = ia+m-1;
%     jalo = ja+i-1;
%     jahi = jalo+ib1-1;
%
%     iblo = ib+i-1;
%     ibhi = ib+m-1;
%     jblo = jb;
%     jbhi = jb+n-1;
%
%     itlo = it+i-1;
%     ithi = itlo+ib1-1;
%     jtlo = jt+i-1;
%     jthi = jtlo+ib1-1;
%
%     ml = m-i+1;
%
%     for ii = 1:k,
%
%        V = tril(A(ialo:iahi,jalo:jahi), -1) + eye(ml,1);
%
%        H = (eye(ml,ml) - V * ( T(itlo:ithi,jtlo:jthi)' * V' ) );
%
%        B(iblo:ibhi,jblo:jbhi) = H*B(iblo:ibhi,jblo:jbhi);
%       
%        ml = ml-1;
%
%        iblo = iblo+1;
%
%        ialo = ialo+1;
%        jalo = jalo+1;
%        jahi = jahi+1;
%        iahi = m;
%
%        itlo = itlo+1;
%        jtlo = jtlo+1;
%        ithi = ithi+1;
%        jthi = jthi+1;
%
%     end

%%%%%%%%%%%%%%%%%%%

     ub = [ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, NaN ];
     ub = [ 17, NaN ];
     ub = [ 10, 7, NaN ];
     ub = [ 9, 2, 4, 2, NaN ];

   tb = [ 4, 3, 13, 2, 4, 8,  7, 9 ];

vb = [ 9, 2, 4, 2 ]


     ll = 1;
     j = 1;
     while( j + tb( ll ) < i )
        j = j + tb( ll );
        ll = ll + 1;
     end
     vb = j + tb( ll ) - i 
     if ( vb > k ) vb = k; end

     ialo = ia+i-1;
     iahi = ia+m-1;
     jalo = ja+i-1;
     jahi = jalo+vb-1;

     iblo = ib+i-1;
     ibhi = ib+m-1;
     jblo = jb;
     jbhi = jb+n-1;

     itlo = it+i-1;
     ithi = itlo+vb-1;
     jtlo = jt+i-1;
     jthi = jtlo+vb-1;

     ml = m-i+1;

     for ii = 1:size(ub,2)-1,

        [jalo:jahi] 

        V = tril(A(ialo:iahi,jalo:jahi), -1) + eye(ml,ub(ii));

        H = (eye(ml,ml) - V * ( T(itlo:ithi,jtlo:jthi)' * V' ) );

        B(iblo:ibhi,jblo:jbhi) = H*B(iblo:ibhi,jblo:jbhi);
       
        ml = ml-ub(ii);

        iblo = iblo+ub(ii);

        ialo = ialo+ub(ii);
        jalo = jalo+ub(ii);
        jahi = jahi+ub(ii+1);
        iahi = m;

        itlo = itlo+ub(ii);
        jtlo = jtlo+ub(ii);
        ithi = ithi+ub(ii+1);
        jthi = jthi+ub(ii+1);

     end

%%%%%%%%%%%%%%%%%%%

end
