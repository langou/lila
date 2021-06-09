%
   close all
   clear
%
%  This is a comparison of recursive DGEQRF and MKL DGEQRF
%
  mt = [ 100, 1000, 5000, 10000, 15000, 20000 ];
  m = 20000;
  n = 20000;
%
%  
%
   mkl_lapack = [ 995.9383, 995.9383, 995.9383,  995.9383,  995.9383,  995.9383 ];
   rec_lapack = [ 998.6091, 924.5512, 991.9947, 1070.4023, 1118.2680, 1194.7758 ];
%   
   mkl_lapack_gflop = [ 21.4203, 21.4203, 21.4203, 21.4203, 21.4203, 21.4203 ];
   rec_lapack_gflop = [ 21.6361, 23.0743, 21.5055, 19.9302, 19.0771, 17.7015 ];
%
%
%
   figure
   ax = axes;
   loglog(mt,mkl_lapack,'LineWidth', 2, 'color', 'b'); hold on;
   loglog(mt,rec_lapack,'LineWidth', 2, 'color', 'r'); hold on;
   grid on;
   axis([ 100 20000 0 1200 ]);
   legend('MKL LAPACK', 'Recursive - LAPACK Panel');
   title({'Profiling of QR Factorization Using LAPACK','Timing Results','Comparing MKL LAPACK to Recursive LAPACK'});
   xlabel('mt - size of block Householder update')
   ylabel('time (seconds)')
   set( ax, 'XTick', mt );
   set( ax, 'XTickLabel', mt );
   set(gca,'FontSize',12);
   grid on
%
%
%
   figure
   ax = axes;
   loglog(mt,mkl_lapack_gflop,'LineWidth', 2, 'color', 'b'); hold on;
   loglog(mt,rec_lapack_gflop,'LineWidth', 2, 'color', 'r'); hold on;
   grid on;
   axis([ 100 20000 0 25 ]);
   legend('MKL LAPACK', 'Recursive - LAPACK Panel');
   title({'Performance Using One Core','m = n = 10,000'});
   xlabel('mt - size of bkock Householder update')
   ylabel('performance for varying mt sizes (GFlop/sec)')
   set( ax, 'XTick', mt );
   set( ax, 'XTickLabel', mt  );
   set(gca,'FontSize',12);
   grid on

   return
