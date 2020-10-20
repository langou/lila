%
   close all
   clear
%
%  This is a comparison of recursive DGEQRF and MKL DGEQRF
%
  mt = [ 10000, 5000, 1000, 100 ];
  m = 10000;
  n = 10000;
%
%  
%
   mkl_lapack = [ 128.0902, 128.0902, 128.0902, 128.0902 ];
   rec_lapack = [ 152.0301, 135.2279, 119.1179, 123.8504 ];
%   
   mkl_lapack_gflop = [ 20.8187, 20.8187, 20.8187, 20.8187 ];
   rec_lapack_gflop = [ 17.5404, 19.7198, 22.3868, 21.5314 ];
%
%
%
   ax = axes;
   loglog(1:4,mkl_lapack,'LineWidth', 2, 'color', 'b'); hold on;
   loglog(1:4,rec_lapack,'LineWidth', 2, 'color', 'r'); hold on;
   grid on;
   axis([ 100 10000 0 153 ]);
   legend('MKL LAPACK', 'Recursive - LAPACK Panel');
   title({'Profiling of QR Factorization Using LAPACK','Time results','Comparing MKL LAPACK to Recursive LAPACK'});
   xlabel('mt - size of block Householder Update')
   ylabel('time (seconds)')
   set( ax, 'XTick', 1:4 );
   set( ax, 'XTickLabel', mt );
   set(gca,'FontSize',12);
   grid on
%
%
%
   figure
   ax = axes;
   plot(1:4,mkl_lapack_gflop,'LineWidth', 2, 'color', 'b'); hold on;
   plot(1:4,rec_lapack_gflop,'LineWidth', 2, 'color', 'r'); hold on;
   grid on;
   axis([ 1 4 0 23 ]);
   legend('MKL LAPACK', 'Recursive - LAPACK Panel');
   title({'Performance Using One Core','m = n = 10,000'});
   xlabel('m')
   ylabel('performance for varying m (GFlop/sec/core)')
   set( ax, 'XTick', 1:4 );
   set( ax, 'XTickLabel', mt  );
   set(gca,'FontSize',12);
   grid on

   return
