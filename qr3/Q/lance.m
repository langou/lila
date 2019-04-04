
unix("gcc -I/Users/dbielich/Documents/lapack.git/CBLAS/include/ -I/Users/dbielich/Documents/lapack.git/LAPACKE/include/ main_flops_xV2N.c flops_xV2N.c flops_syrk.c flops_trmm.c  && source ./test-tmp.sh > K");
load("K");
k=K(:,1);
f=K(:,2);

%f=f -1/6.*(2*k+1).*(k+1).*k + 1/4*k.^2 +5/4*k ;

f = f - (1/3).*(k).*(k).*(k) - (1/3).*(k).*(k)  ;

A = [ k.^3 k.^2 k k.^2.*log2(k) k.*log2(k) log2(k) ones(size(k)) ];
A\f
