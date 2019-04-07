
unix("gcc -I/Users/langou/Desktop/repositories/lapack.git/CBLAS/include/ -I/Users/langou/Desktop/repositories/lapack.git/LAPACKE/include/ main_flops_xV2N.c flops_xV2N.c flops_syrk.c flops_trmm.c  && source ./test-tmp.sh > K");
load("K");
k=K(:,1);
f=K(:,2);

b = f;
b = f - 1/3.*(k).^3 + 1/3.*(k);

[ k f (1/3.*(k).^3 - 1/3.*k) b b./f]

return

%A = [ k k.*ceil(log2(k)) k.*floor(log2(k)) ceil(log2(k)) floor(log2(k)) ones(size(k)) ];
A = [ k.^3 k.^2 k  ones(size(k)) ];
x = A\b

format long g
[ k (b-A*x) ]

