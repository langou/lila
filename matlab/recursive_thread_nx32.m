%
   close all
   clear
%
%  This is a comparison of recursive DGEQRF and MKL DGEQRF
%
  threads = [ 1, 2, 4, 8, 16 ];
  m = 10000;
  n = 10000;
%
%  
%
   mkl_lapack    = [ 128.9239, 65.7663,  35.0425, 19.6928, 11.7264 ];
   rec_lapack10  = [ 220.8271, 120.4531, 66.9560, 51.3747, 40.1148 ];
   rec_lapack100 = [ 124.0594, 66.5030,  35.9387, 20.3060, 13.2543 ];
   rec_lapack500 = [ 118.7819, 62.9670,  33.5687, 18.8811, 11.9339 ];
   rec_lapack1k  = [ 120.3012, 64.1240,  34.1675, 18.9523, 12.2006 ];
   rec_lapack5k  = [ 135.8497, 72.9882,  38.6467, 22.0557, 13.9329 ];
   rec_lapack10k = [ 153.9813, 81.2108,  43.5739, 24.9124, 15.9808 ];

%   
   mkl_lapack_gflop    = [ 20.6840, 40.5476, 76.0981, 135.4135, 227.3878 ];
   rec_lapack_gflop10  = [ 12.0758, 22.1386, 39.8272,  51.9062,  66.4759 ];
   rec_lapack_gflop100 = [ 21.4951, 40.0984, 74.2005, 131.3240, 201.1922 ];
   rec_lapack_gflop500 = [ 22.4501, 42.3502, 79.4391, 141.2351, 223.4538 ];
   rec_lapack_gflop1k  = [ 22.1666, 41.5861, 78.0468, 140.7040, 218.5685 ];
   rec_lapack_gflop5k  = [ 19.6295, 36.5356, 69.0011, 120.9059, 191.3934 ];
   rec_lapack_gflop10k = [ 17.3181, 35.8364, 61.1987, 107.0418, 166.8668 ];
%
%
%
   figure
   ax = axes;
   loglog(threads,mkl_lapack,'LineWidth', 2, 'color', 'k'); hold on;
   loglog(threads,rec_lapack10,'--','LineWidth', 2, 'color', 'r'); hold on;
   loglog(threads,rec_lapack100,'--','LineWidth', 2, 'color', 'b'); hold on;
   loglog(threads,rec_lapack500,'--','LineWidth', 2, 'color', 'm'); hold on;
   loglog(threads,rec_lapack1k,'--','LineWidth', 2, 'color', 'c'); hold on;
   loglog(threads,rec_lapack5k,'--','LineWidth', 2, 'color', 'y'); hold on;
   loglog(threads,rec_lapack10k,'--','LineWidth', 2, 'color', 'g'); hold on;
   grid on;
   axis([ 1 16 0 1000 ]);
   legend('MKL LAPACK', 'mt = 10', 'mt = 100', 'mt = 500', 'mt = 1000', 'mt = 5000', 'mt = 10000');
   title({'Profiling of QR Factorization Using LAPACK','Timing Results','Comparing MKL LAPACK to Recursive LAPACK (nx=32)'});
   xlabel('number of threads')
   ylabel('time (seconds)')
   set( ax, 'XTick', threads );
   set( ax, 'XTickLabel', threads );
   set(gca,'FontSize',12);
   grid on
%
%
%
   figure
   ax = axes;
   plot(threads,mkl_lapack_gflop,'LineWidth', 2, 'color', 'k'); hold on;
   plot(threads,rec_lapack_gflop10,'--','LineWidth', 2, 'color', 'r'); hold on;
   plot(threads,rec_lapack_gflop100,'--','LineWidth', 2, 'color', 'b'); hold on;
   plot(threads,rec_lapack_gflop500,'--','LineWidth', 2, 'color', 'm'); hold on;
   plot(threads,rec_lapack_gflop1k,'--','LineWidth', 2, 'color', 'c'); hold on;
   plot(threads,rec_lapack_gflop5k,'--','LineWidth', 2, 'color', 'y'); hold on;
   plot(threads,rec_lapack_gflop10k,'--','LineWidth', 2, 'color', 'g'); hold on;
   grid on;
   axis([ 1 16 0 230 ]);
   legend('MKL LAPACK', 'mt = 10', 'mt = 100','mt = 500', 'mt = 1000', 'mt = 5000', 'mt = 10000');
   title({'Performance Using One Node','m = n = 10,000, nx = 32, varying mt and number of threads'});
   xlabel('number of threads')
   ylabel('performance for varying mt sizes - (GFlop/sec)')
   set( ax, 'XTick', threads );
   set( ax, 'XTickLabel', threads );
   set(gca,'FontSize',12);
   grid on
%
%
%
   return
