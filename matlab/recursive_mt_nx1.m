%
   close all
   clear
%
%  This is a comparison of recursive DGEQRF and MKL DGEQRF
%
  mt = [ 1, 10, 50, 100:100:1000, 2500, 5000, 10000 ];
  m = 10000;
  n = 10000;
%
%  
%
   mkl_lapack1  = 128.9239*ones(1,size(mt,2));
   rec_lapack1  = [ 978.6006, 220.3600, 134.6966, 124.6622, 119.9866, 119.7910, 118.6702, 119.1469, 118.8736, 118.9859, 119.5035, 120.2647, 119.8432, 125.8140, 136.5313, 153.4843 ];
%
   mkl_lapack2  = 65.7663*ones(1,size(mt,2));
   rec_lapack2  = [ 522.4787, 122.0465, 75.8688, 67.4299, 63.7843, 63.9174, 63.7302, 63.9675, 63.6356, 63.3135, 64.6192, 64.0646, 64.2783, 67.2701, 72.8162, 81.2160 ];
%
   mkl_lapack4  = 35.0425*ones(1,size(mt,2));
   rec_lapack4  = [ 322.3849, 68.0925, 40.3004, 36.3665, 34.3975, 34.3066, 34.3161, 34.0100, 34.2250, 34.0323, 34.9000, 34.3078, 34.6073, 35.9570, 39.1868, 44.0415 ];
%
   mkl_lapack8  = 19.6928*ones(1,size(mt,2));
   rec_lapack8  = [ 276.6118, 52.5733, 22.8147, 20.6116, 19.3580, 19.5373, 19.2606, 19.1295, 19.0925, 19.3535, 19.7195, 19.2968, 19.5098, 20.5141, 21.1619, 25.1728 ];
%
   mkl_lapack16 = 11.7264*ones(1,size(mt,2));
   rec_lapack16 = [ 268.4140, 40.2419, 15.9036, 13.4670, 12.5893, 12.3851, 12.5048, 12.4698, 12.2254, 12.3173, 12.4254, 12.3331, 12.4106, 13.0019, 14.3223, 16.5065 ];
%   
   mkl_lapack_gflop1  = 20.6840*ones(1,size(mt,2));
   rec_lapack_gflop1  = [ 2.7250, 12.1014, 19.7976, 21.3911, 22.2247, 22.2610, 22.4712, 22.3813, 22.4238, 22.4116, 22.3145, 22.1733, 22.2513, 21.1953, 19.5315, 17.3743 ];
%   
   mkl_lapack_gflop2  = 40.5476*ones(1,size(mt,2));
   rec_lapack_gflop2  = [ 5.1039, 21.8496, 35.6178, 39.5472, 41.8076, 41.7205, 41.8431, 41.6879, 41.9053, 42.1185, 41.2674, 41.6247, 41.4863, 39.6412, 36.6219, 32.8343 ];
%   
   mkl_lapack_gflop4  = 76.0981*ones(1,size(mt,2));
   rec_lapack_gflop4  = [ 8.2717, 39.1624, 66.1697, 73.3275, 77.5251, 77.7305, 77.7089, 78.4084, 77.9159, 78.3570, 76.4088, 77.7277, 77.0551, 74.1628, 68.0501, 60.5490 ];
%   
   mkl_lapack_gflop8  = 135.4135*ones(1,size(mt,2));
   rec_lapack_gflop8  = [ 9.6405, 50.7229, 116.8839, 129.3769, 137.7552, 136.4914, 138.4521, 139.4006, 139.6707, 137.7873, 135.2299, 138.1920, 136.6838, 129.9918, 121.4227, 105.9344 ];
%   
   mkl_lapack_gflop16 = 227.3878*ones(1,size(mt,2));
   rec_lapack_gflop16 = [ 9.9349, 66.2659, 167.6767, 198.0147, 211.8197, 215.3118, 213.2513, 213.8505, 218.1250, 216.4977, 214.6140, 216.2196, 214.8609, 205.0989, 186.1893, 161.5530 ];
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
   title({'Performance Using One Node','m = n = 10,000'});
   xlabel('mt - size of block Householder update')
   ylabel('performance for varying mt sizes (GFlop/sec)')
   set( ax, 'XTick', mt );
   set( ax, 'XTickLabel', mt  );
   set(gca,'FontSize',12);
   grid on

   return
