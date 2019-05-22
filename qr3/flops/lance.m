clear
unix("make  && source ./test-flops1.sh > K");
load("K");
m=K(:,1);
n=K(:,2);
f=K(:,3);
b = 1/3 * m.*n.^2 - m.*n + 2/3* m -2/21*n.^3  +1/2 *n.^2 -5/6* n +3/7  ;

[ f b f-b]

A = [ m.*n.^2 m.*n m n.^3 n.^2 n ones(size(n)) ]


A\(f-b)
