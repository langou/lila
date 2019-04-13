
unix("make  && source ./test-tmp.sh > K");
load("K");
m=K(:,1);
n=K(:,2);
f=K(:,4);

b = f;
k = m.^2.*n + m.*n;

[ b k b-k ]

%A = [ m.^2.*n + m.*n  ];
%format rat
%x = A\b

%format long g
%[ k b round(b) A*x b-A*x ]

