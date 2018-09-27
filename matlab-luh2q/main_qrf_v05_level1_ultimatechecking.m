%
   clear
%
   fprintf('\n');
%
   threshold = 1e-12;
   extrarows_max = 1;
   mt_max = 50;
   log10KA_max = 1.50;
   nb_blocks_lvl1_max = 5;
   nb_blocks_lvl2_max = 6;
   sz_blocks_lvl3_max = 10;
   nbtest = 10;
%
   fprintf('        m    n    mt \n');
%
   for test = 1:nbtest,
   clear global
%
   global nb_lvl1;
   global ii_lvl2;
   global nb_lvl2;
%
   mt = ceil(rand(1)*mt_max);
   nb_blocks_lvl1 = ceil(rand(1)*nb_blocks_lvl1_max);
   log10KA = rand(1)*log10KA_max;
   nb_lvl2 = cell(1,nb_blocks_lvl1);
   for i = 1:nb_blocks_lvl1, 
   nb_blocks_lvl2 = ceil(rand(1)*nb_blocks_lvl2_max);
   for j = 1:nb_blocks_lvl2, 
   nb_lvl2{i}(1,j) = ceil(rand(1)*sz_blocks_lvl3_max);
   end
   end
%
   nb_lvl1 = zeros(1,nb_blocks_lvl1);
   nb_blocks_lvl2 = zeros(1,nb_blocks_lvl1);
   for i = 1:nb_blocks_lvl1,
      nb_lvl1(i) = sum( nb_lvl2{i} );
      nb_blocks_lvl2(i) = size(nb_lvl2{i},2);
   end
   n = sum(nb_lvl1);
%
   m = n + ( ceil(rand(1)*(extrarows_max+1))-1 );
%
   fprintf(' %4d %4d %4d %4d',test,m,n,mt);
%
   [check] = main_qrf_v05_level1_ultimatechecking_do( m, n, mt, log10KA );
%
   fprintf(' %6.1e', check(1));
   fprintf(' %6.1e', check(2));
   fprintf(' %6.1e', check(3));
%  fprintf(' %6.1e', check(4)); % this check will be off for square matrices for some reason we understand
   fprintf(' %6.1e', check(5));
   fprintf(' %6.1e', check(6));
%
   for i = 1:nb_blocks_lvl1, 
   for j = 1:size(nb_lvl2{i},2), 
   fprintf(' %4d ',i );
   end
   end
   fprintf('\n');
%
   fprintf('                                                            ');
%
   for i = 1:nb_blocks_lvl1, 
   for j = 1:size(nb_lvl2{i},2), 
   fprintf(' %4d ',nb_lvl2{i}(j) );
   end
   end
   fprintf('\n');
%
%  note that we do not check against #4
   if(( check(1) > threshold ) || ( check(2) > threshold ) || ( check(3) > threshold ) || ( check(5) > threshold ) || ( check(6) > threshold ) ), fprintf('       ^^^   shoooooooo  ^^^\n'); return; end
%
   end
