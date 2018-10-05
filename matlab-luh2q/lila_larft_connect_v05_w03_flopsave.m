   function [ T ] = lila_larft_connect_v05_w03_flopsave( m, n, i, mt, A, T )
%
   vb= mt - mod(i-1,mt) ;
   if ( vb > n ), vb = n; end
%
   itlo = mod(i-1,mt)+1;
   ithi = itlo+vb-1;
   jtlo = i;
   jthi = jtlo+vb-1;
%
%   T(1:itlo-1,jtlo:jthi) = ...
%       - ( T(1:itlo-1,jtlo-(itlo)+1:jtlo-1) * A(jtlo:m,jtlo-(itlo)+1:jtlo-1)' )* ...
%           ( (tril(A(jtlo:m,jtlo:jthi),-1) + eye(size((A(jtlo:m,jtlo:jthi)))) ) * T(itlo:ithi,jtlo:jthi) );
%
%  work = ( A(jtlo:m,jtlo-(itlo)+1:jtlo-1)' ) * (tril(A(jtlo:m,jtlo:jthi),-1) + eye(size((A(jtlo:m,jtlo:jthi)))) );    
%  T(1:itlo-1,jtlo:jthi) = - ( T(1:itlo-1,jtlo-(itlo)+1:jtlo-1) * work * T(itlo:ithi,jtlo:jthi) );
%
%   T(1:itlo-1,jtlo:jthi) = ( A(jtlo:jthi,jtlo-itlo+1:jtlo-1)' ) * (tril(A(jtlo:jthi,jtlo:jthi),-1) + eye(vb) );    
%   T(1:itlo-1,jtlo:jthi) = T(1:itlo-1,jtlo:jthi) + ( A(jthi+1:m,jtlo-(itlo)+1:jtlo-1)' ) * A(jthi+1:m,jtlo:jthi);    
%   T(1:itlo-1,jtlo:jthi) = - T(1:itlo-1,jtlo-(itlo)+1:jtlo-1) * T(1:itlo-1,jtlo:jthi) ;
% T(1:itlo-1,jtlo:jthi) = T(1:itlo-1,jtlo:jthi) * T(itlo:ithi,jtlo:jthi);
%
    lda = -1;
    ldt = -1;
%
   T(1:itlo-1,jtlo:jthi) = A(jtlo:jthi,jtlo-itlo+1:jtlo-1)';
   [ T ] = blas_trmm ( 'R', 'L', 'N', 'U', itlo-1, vb, (+1.0e+00), A, jtlo, jtlo-itlo+1, ldt, T, 1, jtlo, ldt );
   [ T ] = blas_gemm ( 'T', 'N', itlo-1, vb, m-jthi, (+1.0e+00), A, jthi+1, jtlo-(itlo)+1, lda, A, jthi+1, jtlo, lda, (+1.0e+00), T, 1, jtlo, ldt );
%%%%%%
   [ T ] = blas_trmm ( 'L', 'U', 'N', 'N', itlo-1, vb, (-1.0e+00), T, 1, jtlo-(itlo)+1, ldt, T, 1, jtlo, ldt );
%  T(1:itlo-1,jtlo:jthi) = - T(1:itlo-1,jtlo-(itlo)+1:jtlo-1) * T(1:itlo-1,jtlo:jthi) ;
%
   [ T ] = blas_trmm ( 'R', 'U', 'N', 'N', itlo-1, vb, (+1.0e+00), T, itlo, jtlo, ldt, T, 1, jtlo, ldt );
%     T(1:itlo-1,jtlo:jthi) = T(1:itlo-1,jtlo:jthi) * T(itlo:ithi,jtlo:jthi);
%%%%%%



   end
