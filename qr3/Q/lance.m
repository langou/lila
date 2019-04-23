clear
unix("make  && source ./test-tmp4.sh > K");
load("K");
m=K(:,1);
n=K(:,2);
k=K(:,3);

f= K(:,7) ;

format long g
A = [ m.*n m.*k k.*n m n k  ones(size(n,1))];

x = A\f

[ m n k f A*x f-A*x ]

