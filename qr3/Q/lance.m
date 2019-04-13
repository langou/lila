
unix("make  && source ./test-tmp.sh > K");
load("K");
n=K(:,1);
f=K(:,3);

b = f;
k = [ (2*n.^3+3*n.^2-5*n)/6   ];

[ n b k b-k ]
return

A = [ n.^3 n.^2 n ones(size(n))  ];
format rat
x = A\b


format long g
[ k b round(b) A*x b-A*x ]

