clear
unix("make  && source ./test.sh > K");
load("K");
m=K(:,1);
n=K(:,2);
b=K(:,3);

f= K(:,6) ;

format long g

nb = b;
kb = floor( n ./ b );

A= [
 n.*nb.*kb, ...
 n.*nb, ...
 n.*kb, ...
 nb.*nb.*kb.*kb, ...
 nb.*kb, ...
 nb.*nb, ...
 kb.*kb, ...
 nb.*nb.*kb, ...
 nb.*kb.*kb, ...
 ones(size(n)) ];

x = A\f

[ m n nb f A*x f-A*x ];

