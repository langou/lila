%
   clear
%
   m = 20;
   n1 = 7;
   n2 = 5;
%
   n = n1+n2;
   if ( m < n ) fprintf('m < n\n'); return; end;
%
   A = randn(m,n);
%
   A1 = A(1:m,1:n1);
   A2 = A(1:m,n1+1:n);
%
%%%%%%%%

   [ Q1, R11, V1, T11 ] = geqr3wQ ( A1 );

   [ A2 ] = larfb( V1, T11, A2 );
   R12 = A2(1:n1,1:n2);

   [ Q2, R22, V2, T22 ] = geqr3wQ ( A2(n1+1:m,1:n2) );
   V2 = [ zeros(n1,n2); V2 ];

oo = 5;
if ( oo == 0 )

   Q2 = [ zeros(n1,n2); Q2 ];
   Q2 = Q2 - V1 * ( T11 * ( V1' * Q2 ) ) ;
   T12 = - ( T11 * ( V1(n1+1:m,1:n1)' * V2(n1+1:m,1:n2) ) * T22 ) ;

elseif ( oo == 1 )

   Q2 = [ zeros(n1,n2); eye(n2,n2); zeros(m-n1-n2,n2)];
   Q2 = Q2 - V2 * ( T22 * ( V2' * Q2 ) ) ;
   Q2 = Q2 - V1 * ( T11 * ( V1' * Q2 ) ) ;
   T12 = - ( T11 * ( V1(n1+1:m,1:n1)' * V2(n1+1:m,1:n2) ) * T22 ) ;

elseif ( oo == 2 )

   Q2 = [ zeros(n1,n2); eye(n2,n2); zeros(m-n1-n2,n2)];
   Q2 = ( Q2 - V2 * ( T22 * ( V2' * Q2 ) ) ) - V1 * ( T11 * ( V1' * ( Q2 - V2 * ( T22 * ( V2' * Q2 ) ) ) ) ) ;
   T12 = - ( T11 * ( V1(n1+1:m,1:n1)' * V2(n1+1:m,1:n2) ) * T22 ) ;

elseif ( oo == 3 )

   Q2 = [ zeros(n1,n2); eye(n2,n2); zeros(m-n1-n2,n2)];
   Q12 = [ zeros(n1,n2) ];
   Q22 = [ eye(n2,n2) ];
   Q32 = [ zeros(m-n1-n2,n2) ];

   V11 = [ V1(1:n1,1:n1) ];
   V21 = [ V1(n1+1:n,1:n1) ];
   V31 = [ V1(n+1:m,1:n1) ];

   V12 = [ V2(1:n1,1:n2) ];
   V22 = [ V2(n1+1:n,1:n2) ];
   V32 = [ V2(n+1:m,1:n2) ];

   Q12 = ( Q12 - V12 * ( T22 * ( V2' * Q2 ) ) ) - V11 * ( T11 * ( V1' * ( Q2 - V2 * ( T22 * ( V2' * Q2 ) ) ) ) ) ;
   Q22 = ( Q22 - V22 * ( T22 * ( V2' * Q2 ) ) ) - V21 * ( T11 * ( V1' * ( Q2 - V2 * ( T22 * ( V2' * Q2 ) ) ) ) ) ;
   Q32 = ( Q32 - V32 * ( T22 * ( V2' * Q2 ) ) ) - V31 * ( T11 * ( V1' * ( Q2 - V2 * ( T22 * ( V2' * Q2 ) ) ) ) ) ;

   Q2 = [ Q12; Q22; Q32 ];

   T12 = - ( T11 * ( V21' * V22 + V31'*V32 ) * T22 ) ;

elseif ( oo == 4 )

   V11 = [ V1(1:n1,1:n1) ];
   V21 = [ V1(n1+1:n,1:n1) ];
   V31 = [ V1(n+1:m,1:n1) ];

   V22 = [ V2(n1+1:n,1:n2) ];
   V32 = [ V2(n+1:m,1:n2) ];

   Q12 =  - V11 * T11 * (  V21' - ( ( V21' * V22 + V31' * V32  ) *  T22 * V22'  ) ) ;
   Q22 = eye(n2,n2) - V22 * T22 * V22' - V21 * ( T11 * (  V21' - ( V21' * V22 ) *  T22 * V22' - ( V31' * V32 ) * T22 * V22' ) ) ;
   Q32 = ( - V32 * T22 * V22' ) - V31 * ( T11 * (  V21' - ( V21' * V22 ) *  T22 * V22' - ( V31' * V32 ) * T22 * V22' ) ) ;

   Q2 = [ Q12; Q22; Q32 ];

   T12 = - ( T11 * ( V21' * V22 + V31' * V32 ) * T22 ) ;

elseif ( oo == 5 )

   V11 = [ V1(1:n1,1:n1) ];
   V21 = [ V1(n1+1:n,1:n1) ];
   V31 = [ V1(n+1:m,1:n1) ];

   V22 = [ V2(n1+1:n,1:n2) ];
   V32 = [ V2(n+1:m,1:n2) ];

   T12 = - ( T11 * ( V21' * V22 + V31' * V32 ) * T22 ) ;

   Q12 = - ( T11 *  V21' + T12 * V22' ) ;

   Q22 = - T22 * V22' ;

   Q32 = V32 * Q22 + V31 * Q12 ;

   Q22 = eye(n2,n2) + V22 * Q22 + V21 * Q12 ;

   Q12 = V11  * Q12 ;

   Q2 = [ Q12; Q22; Q32 ];

end

   Q = [Q1 Q2];
   V = [V1 V2];
   R = [ R11 R12; zeros(n2,n1) R22 ];
   T = [ T11 T12; zeros(n2,n1) T22 ];

%  [ Q, R, V, T ] = geqr3wQ ( A );

   H = ( eye(m) - V * T * V' );
   fprintf('|| A - Q*R || / || A ||          = %6.1e\n',norm(A-Q*R,'fro')/norm(A,'fro'));
   fprintf('|| I - Q''*Q ||                   = %6.1e\n',norm(eye(n) - Q'*Q,'fro'));
   fprintf('|| tril( H''*A ) || / || A ||     = %6.1e\n',norm(tril(H'*A,-1),'fro')/norm(A,'fro'));
   fprintf('|| triu( H''*A ) - R || / || A || = %6.1e\n',norm([eye(n);zeros(m-n,n)]'*triu(H'*A,-1)-R,'fro')/norm(A,'fro'));
   fprintf('|| H''*Q - I ||                   = %6.1e\n',norm(H'*Q-eye(m,n),'fro'));


