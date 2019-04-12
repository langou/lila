
unix("make  && source ./test-tmp.sh > K");
load("K");
k=K(:,1);
f=K(:,2);

b = f;
b = f - 2/3.*(k).^3 + 2/3 * k;


%[ k f (2/3.*(k).^3 - 2/3.*k) b b./f]

%return

%A = [ k k.*ceil(log2(k)) k.*floor(log2(k)) ceil(log2(k)) floor(log2(k)) ones(size(k)) ];
A = [ k.^3 k.^2 k  ones(size(k)) ];
%format rat
x = A\b

format long g
[ k b round(b) A*x b-A*x ]

