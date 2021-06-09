%
   close all
   clear
%
%  This is a comparison of recursive DGEQRF and MKL DGEQRF
%
  mt = [ 10, 100, 500, 1000, 5000, 10000 ];
  m = 10000;
  n = 10000;
%
%  
%
   mkl_lapack1  = 128.9239*ones(1,size(mt,2));
   rec_lapack1  = [ 220.8271, 124.0594, 118.7819, 120.3012, 135.8497, 153.9813 ];
%
   mkl_lapack2  = 65.7663*ones(1,size(mt,2));
   rec_lapack2  = [ 120.4531, 66.5030, 62.9670, 64.1240, 72.9882, 81.2108 ];
%
   mkl_lapack4  = 35.0425*ones(1,size(mt,2));
   rec_lapack4  = [ 66.9560, 35.9387, 33.5687, 34.1675, 38.6467, 43.5739 ];
%
   mkl_lapack8  = 19.6928*ones(1,size(mt,2));
   rec_lapack8  = [ 51.3747, 20.3060, 18.8811, 18.9523, 22.0557, 24.9124 ];
%
   mkl_lapack16 = 11.7264*ones(1,size(mt,2));
   rec_lapack16 = [ 40.1148, 13.2543, 11.9339, 12.2006, 15.9808, 15.9808 ];
%   
   mkl_lapack_gflop1  = 20.6840*ones(1,size(mt,2));
   rec_lapack_gflop1  = [ 12.0758, 21.4951, 22.4501, 22.1666, 19.6295, 17.3181 ];
%   
   mkl_lapack_gflop2  = 40.5476*ones(1,size(mt,2));
   rec_lapack_gflop2  = [ 22.1386, 40.0984, 42.3502, 41.5861, 36.5356, 32.8364 ];
%   
   mkl_lapack_gflop4  = 76.0981*ones(1,size(mt,2));
   rec_lapack_gflop4  = [ 39.8272, 74.2005, 79.4391, 78.0468, 69.0011, 61.1987 ];
%   
   mkl_lapack_gflop8  = 135.4135*ones(1,size(mt,2));
   rec_lapack_gflop8  = [ 51.9062, 131.3240, 141.2351, 140.7040, 120.9059, 107.0418 ];
%   
   mkl_lapack_gflop16 = 227.3878*ones(1,size(mt,2));
   rec_lapack_gflop16 = [ 66.4759, 201.1922, 223.4538, 218.5685, 191.3934, 166.8668 ];
%
%
%
   figure
   ax = axes;
   loglog(mt,mkl_lapack1,'LineWidth', 2, 'color', 'b'); hold on;
   loglog(mt,rec_lapack1,'--','LineWidth', 2, 'color', 'b'); hold on;
   loglog(mt,mkl_lapack2,'LineWidth', 2, 'color', 'r'); hold on;
   loglog(mt,rec_lapack2,'--','LineWidth', 2, 'color', 'r'); hold on;
   loglog(mt,mkl_lapack4,'LineWidth', 2, 'color', 'c'); hold on;
   loglog(mt,rec_lapack4,'--','LineWidth', 2, 'color', 'c'); hold on;
   loglog(mt,mkl_lapack8,'LineWidth', 2, 'color', 'g'); hold on;
   loglog(mt,rec_lapack8,'--','LineWidth', 2, 'color', 'g'); hold on;
   loglog(mt,mkl_lapack16,'LineWidth', 2, 'color', 'k'); hold on;
   loglog(mt,rec_lapack16,'--','LineWidth', 2, 'color', 'k'); hold on;
   grid on;
   axis([ 1 10000 0 1000 ]);
   legend('MKL LAPACK - 1 thread', 'Recursive - LAPACK Panel - 1 thread', 'MKL LAPACK - 2 threads', 'Recursive - LAPACK Panel - 2 threads', 'MKL LAPACK - 4 threads', 'Recursive - LAPACK Panel - 4 threads', 'MKL LAPACK - 8 threads', 'Recursive - LAPACK Panel - 8 threads', 'MKL LAPACK - 16 threads', 'Recursive - LAPACK Panel - 16 threads');
   title({'Profiling of QR Factorization Using LAPACK','Timing Results','Comparing MKL LAPACK to Recursive LAPACK (nx=32)'});
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
   plot(mt,mkl_lapack_gflop1,'LineWidth', 2, 'color', 'b'); hold on;
   plot(mt,rec_lapack_gflop1,'--','LineWidth', 2, 'color', 'b'); hold on;
   plot(mt,mkl_lapack_gflop2,'LineWidth', 2, 'color', 'r'); hold on;
   plot(mt,rec_lapack_gflop2,'--','LineWidth', 2, 'color', 'r'); hold on;
   plot(mt,mkl_lapack_gflop4,'LineWidth', 2, 'color', 'c'); hold on;
   plot(mt,rec_lapack_gflop4,'--','LineWidth', 2, 'color', 'c'); hold on;
   plot(mt,mkl_lapack_gflop8,'LineWidth', 2, 'color', 'g'); hold on;
   plot(mt,rec_lapack_gflop8,'--','LineWidth', 2, 'color', 'g'); hold on;
   plot(mt,mkl_lapack_gflop16,'LineWidth', 2, 'color', 'k'); hold on;
   plot(mt,rec_lapack_gflop16,'--','LineWidth', 2, 'color', 'k'); hold on;
   grid on;
   axis([ 1 10000 0 230 ]);
   legend('MKL LAPACK - 1 thread', 'Recursive - LAPACK Panel - 1 thread', 'MKL LAPACK - 2 threads', 'Recursive - LAPACK Panel - 2 threads', 'MKL LAPACK - 4 threads', 'Recursive - LAPACK Panel - 4 threads', 'MKL LAPACK - 8 threads', 'Recursive - LAPACK Panel - 8 threads', 'MKL LAPACK - 16 threads', 'Recursive - LAPACK Panel - 16 threads');
   title({'Performance Using One Node','m = n = 10,000, nx = 32'});
   xlabel('mt - size of block Householder update')
   ylabel('performance for varying mt sizes (GFlop/sec)')
   set( ax, 'XTick', mt );
   set( ax, 'XTickLabel', mt  );
   set(gca,'FontSize',12);
   grid on

   return
