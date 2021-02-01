%
   close all
   clear
%
%  This is a comparison of recursive DGEQRF and MKL DGEQRF
%
  mt = [ 100, 1000, 5000, 10000, 15000 ];
  m = 15000;
  n = 15000;
%
%  
%
   mkl_lapack = [ 421.7878, 421.7878, 421.7878, 421.7878, 421.7878 ];
   rec_lapack = [ 420.3326, 393.9459, 426.1572, 461.6981, 506.8162 ];
%   
   mkl_lapack_gflop = [ 21.3377, 21.3377, 21.3377, 21.3377, 21.3377 ];
   rec_lapack_gflop = [ 21.4116, 22.8458, 21.1190, 19.4933, 17.7580 ];
%
%
%
   figure
   ax = axes;
   loglog(mt,mkl_lapack,'LineWidth', 2, 'color', 'b'); hold on;
   loglog(mt,rec_lapack,'LineWidth', 2, 'color', 'r'); hold on;
   grid on;
   axis([ 100 15000 0 510 ]);
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
   axis([ 100 15000 0 25 ]);
   legend('MKL LAPACK', 'Recursive - LAPACK Panel');
   title({'Performance Using One Core','m = n = 10,000'});
   xlabel('mt - size of bkock Householder update')
   ylabel('performance for varying mt sizes (GFlop/sec)')
   set( ax, 'XTick', mt );
   set( ax, 'XTickLabel', mt  );
   set(gca,'FontSize',12);
   grid on

   return
