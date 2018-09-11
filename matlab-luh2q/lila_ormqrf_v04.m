   function [ B ] = lila_ormqrf_v04( m, n, k, i, A, ia, ja, lda, B, ib, jb, ldb, T, it, jt, ldt, tb )

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

%    ub = [ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, NaN ];
%    ub = [ 17, NaN ];
%    ub = [ 10, 7, NaN ];
%    ub = [ 9, 2, 4, 2, NaN ];

%  tb = [ 4, 3, 13, 2, 4, 8,  7, 9 ];
%
%vb = [ 9, 2, 4, 2 ]
%

     ll = 1;
     jj = 1;
     while( jj + tb( ll ) < i )
        jj = jj + tb( ll );
        ll = ll + 1;
     end
     vb = jj + tb( ll ) - i;
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

%ii = 0;

     not_done = 1;

     jj = 1;

%    while ( ii < 4 )

     while ( not_done == 1 )

%       [jalo:jahi] 

        V = tril(A(ialo:iahi,jalo:jahi), -1) + eye(ml,vb);

        H = (eye(ml,ml) - V * ( T(itlo:ithi,jtlo:jthi)' * V' ) );

        B(iblo:ibhi,jblo:jbhi) = H*B(iblo:ibhi,jblo:jbhi);
       
        if ( jj + vb - 1 == k ) 

           not_done = 0; 

        else if ( jj + vb - 1 > k ) 

            fprintf('defensive programming, we should never have been there, abort\n'); return;

        else

           ml = ml-vb;

           ialo = ialo+vb;
           iblo = iblo+vb;
           jalo = jalo+vb;
           itlo = itlo+vb;
           jtlo = jtlo+vb;

           jj = jj + vb;

           ll = ll + 1;

           if ( ( jj+tb(ll)-1 ) <= k ), vb = tb(ll); else vb = k-jj+1; end;

           jahi = jahi+vb;
           ithi = ithi+vb;
           jthi = jthi+vb;


end


     end

%%%%%%%%%%%%%%%%%%%

end
