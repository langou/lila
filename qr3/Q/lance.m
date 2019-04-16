
unix("make  && source ./test-tmp.sh > K");
load("K");
m=K(:,1);
n=K(:,2);
f=K(:,6);

b=          2.0000*n.^2.*m    +2*n.*m -2/3 *n.^3       -4/3*n     ;

%format long g
[ m n b f b-f ]
return


%A = [   n.^3 n  n.^2.*m n.*m   ];

A = [ m.^3 m.^2 m n.^3 n.^2 n m.^2.*n n.^2.*m n.*m ones(size(n))  ];


  %m.^3		m.^2 	m	 n.^3		 n.^2 n   k.^3      k.^2      k         m.^2.*n   m.^2.*k   n.^2.*m    n.^2.*k  k.^2.*m    k.^2.*n   n.*m      k.*m      k.*n      ones(size(n))  ];
  %-0.0000    	0.0000  -0.0000   -0.6667   	 0.0000   -1.3333   -0.0000   -0.0000   -0.0000   -0.0000   -0.0000    2.0000   -0.0000    0.0000    0.0000    2.0000    0.0000    0.0000    .
%format rat
x = A\f


%format long g
[ m n f round(f) A*x f-A*x ]

