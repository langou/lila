clear
unix("make  && source ./test-flops1.sh > K");
load("K");
m=K(:,1);
n=K(:,2);
f=K(:,3);
%b = 1/3 * m.*n.^2 - m.*n + 2/3* m -2/21*n.^3  +1/2 *n.^2 -5/6* n +3/7  ;
b = zeros(size(f));
%b = 7/3 * m.*n.^2  + 2/3* m -9/14*n.^3   +11/2* n +1/7  ;
b = ( 98 * m.*n.^2  + 28* m -27*n.^3   +231* n +6 ) / 42  ;

[ f b f-b]
return

A = [ m.*n.^2 m.*n m n.^3 n.^2 n ones(size(n)) ]


A\(f-b)
