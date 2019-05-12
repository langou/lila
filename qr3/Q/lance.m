clear
unix("make  && source ./test-flops.sh > K");
load("K");
m=K(:,1);
n=K(:,2);
f=K(:,3);

%g = 1/6 * n.^3 - 1/2 * n.^2 + 1/3 .* n;         % 2^k
%g = 1/6 * n.^3 - 2/3 * n.^2 + 1/6 .* n + 1;     % 2^k+1
%g = 1/6 * n.^3 - 7/12 * n.^2 + 1/6 .* n + 2;    % 2^k+2
%g = 1/6 * n.^3 - 11/16 * n.^2 + 5/24 .* n + 49/16;  % 2^k+3
%g = 1/6 * n.^3 - 13/24 * n.^2 + 1/6 .* n + 4;  % 2^k+4

Y = [ m n (f-1/6 * n.^3) ]

A = [ ones(size(n)) n  n.^(3/2) n.^2 n.^(5/2) ];

x = A\(f-1/6 * n.^3)

 flops = flops + 2^(k-1) * (16-1)*16*16 + 2 * (8-1)*8*8 + 4 * (4-1)*4*4  + 2^(k-1) * (2-1)*2*2
flops = 
for k = 1:log2(32)-1,
z = n/2;
2^(k-1)
 flops = flops + 2^(k-1) * (16-1)*16*16 + 2 * (8-1)*8*8 + 4 * (4-1)*4*4  + 2^(k-1) * (2-1)*2*2
end
 1 1 0
 2 2 4
 1 1 0
 4 4 48
 1 1 0
 2 2 4
 1 1 0
 8 8 448
 1 1 0
 2 2 4
 1 1 0
 4 4 48
 1 1 0
 2 2 4
 1 1 0
 16 16 3840
 1 1 0
 2 2 4
 1 1 0
 4 4 48
 1 1 0
 2 2 4
 1 1 0
 8 8 448
 1 1 0
 2 2 4
 1 1 0
 4 4 48
 1 1 0
 2 2 4
 1 1 0


return

 

bin = zeros(size(n,1),4);
for i=1:size(n,1),
	ii = i;
%	if ( ii >= 32 ) ii = ii - 32; bin(i,1) = 32;  end;
%	if ( ii >= 16 ) ii = ii - 16; bin(i,2) = 16;  end;
%	if ( ii >=  8 ) ii = ii -  8; bin(i,3) =  8;  end;
%	if ( ii >=  4 ) ii = ii -  4; bin(i,4) =  4;  end;
%	if ( ii >=  2 ) ii = ii -  2; bin(i,5) =  2;  end;
%	if ( ii >=  1 ) ii = ii -  1; bin(i,6) =  1;  end;

	if ( ii >= 32 ) ii = ii - 32; bin(i,1) =  1;  end;
	if ( ii >= 16 ) ii = ii - 16; bin(i,2) =  2;  end;
	if ( ii >=  8 ) ii = ii -  8; bin(i,3) =  4;  end;
	if ( ii >=  4 ) ii = ii -  4; bin(i,4) =  8;  end;
	if ( ii >=  2 ) ii = ii -  2; bin(i,5) = 16;  end;
	if ( ii >=  1 ) ii = ii -  1; bin(i,6) = 32;  end;

end
return

z = f-g;

A = [ ones(size(n)) bin ];

x = A\(f)

bin(30:end-1,:) \ f(30:end-1,1)

% 54 = 32 + 16 + 4 + 2 => 270
% 56 = 32 + 16 + 8 => 56
% 48 = 32 + 16 => 16
% 40 = 32 +  8 => 32
% 36 = 32 +  4 => 56
% 34 = 32 +  2 => 100

% 32 * 0 + 16 *1 + 8 * 4

% 32 *0 + 2 * 33



