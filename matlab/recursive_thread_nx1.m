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
   mkl_lapack    = [ 128.9239,  65.7663,  35.0425,  19.6928,  11.7264 ];
   rec_lapack1   = [ 978.6006, 522.4787, 322.3849, 276.6118, 268.3149 ];
   rec_lapack10  = [ 220.3600, 122.0465,  68.0925,  52.5733,  40.2419 ];
   rec_lapack50  = [ 134.6966,  74.8688,  40.3004,  22.8147,  15.9036 ];
   rec_lapack100 = [ 124.6622,  67.4299,  36.3665,  20.6116,  13.4670 ];
   rec_lapack200 = [ 119.9866,  63.7843,  34.3975,  19.3580,  12.5893 ];
   rec_lapack300 = [ 119.7910,  63.9174,  34.3066,  19.5373,  12.3851 ];
   rec_lapack400 = [ 118.6702,  63.7302,  34.3161,  19.2606,  12.5048 ];
   rec_lapack500 = [ 119.1469,  63.9675,  34.0100,  19.1295,  12.4698 ];
   rec_lapack600 = [ 118.8736,  63.6356,  34.2250,  19.0925,  12.2254 ];
   rec_lapack700 = [ 118.9859,  63.3135,  34.0323,  19.3535,  12.3173 ];
   rec_lapack800 = [ 119.5035,  64.6192,  34.9000,  19.7195,  12.4254 ];
   rec_lapack900 = [ 120.2647,  64.0646,  34.3078,  19.2968,  12.3331 ];
   rec_lapack1k  = [ 119.8432,  64.2783,  34.6073,  19.5098,  12.4106 ];
   rec_lapack2k  = [ 125.8140,  67.2701,  35.9570,  20.5141,  13.0019 ];
   rec_lapack5k  = [ 136.5313,  72.8162,  39.1868,  21.9619,  14.3223 ];
   rec_lapack10k = [ 153.4843,  81.2160,  44.0415,  25.1728,  16.5065 ];

%   
   mkl_lapack_gflop    = [ 20.6840, 40.5476, 76.0981, 135.4135, 227.3878 ];
   rec_lapack_gflop1   = [  2.7250,  5.1039,  8.2717,   9.6405,   9.9349 ];
   rec_lapack_gflop10  = [ 12.1014, 21.8496, 39.1624,  50.7229,  66.2659 ];
   rec_lapack_gflop50  = [ 19.7976, 35.6178, 66.1697, 116.8839, 167.6767 ];
   rec_lapack_gflop100 = [ 21.3911, 39.5472, 73.3275, 129.3769, 198.0147 ];
   rec_lapack_gflop200 = [ 22.2247, 41.8076, 77.5251, 137.7552, 211.8197 ];
   rec_lapack_gflop300 = [ 22.2610, 41.7205, 77.7305, 136.4914, 215.3118 ];
   rec_lapack_gflop400 = [ 22.4712, 41.8431, 77.7089, 138.4521, 213.2513 ];
   rec_lapack_gflop500 = [ 22.3813, 41.6879, 78.4084, 139.4006, 213.8505 ];
   rec_lapack_gflop600 = [ 22.4238, 41.9053, 77.9159, 139.6707, 218.1250 ];
   rec_lapack_gflop700 = [ 22.4116, 42.1185, 78.3570, 137.7873, 216.4977 ];
   rec_lapack_gflop800 = [ 22.3145, 41.2674, 76.4088, 135.2299, 214.6140 ];
   rec_lapack_gflop900 = [ 22.1733, 41.6247, 77.7277, 138.1920, 216.2196 ];
   rec_lapack_gflop1k  = [ 22.2513, 41.4863, 77.0551, 136.6838, 214.3878 ];
   rec_lapack_gflop2k  = [ 21.1953, 39.6412, 74.1628, 129.9918, 205.0989 ];
   rec_lapack_gflop5k  = [ 19.5315, 36.6219, 68.0501, 121.4227, 186.1893 ];
   rec_lapack_gflop10k = [ 17.3743, 32.8343, 60.5490, 105.9344, 161.5530 ];
%
%
%
   figure
   ax = axes;
   loglog(threads,mkl_lapack,'LineWidth', 2, 'color', 'k'); hold on;
   loglog(threads,rec_lapack1,'--','LineWidth', 2, 'color', 'r'); hold on;
   loglog(threads,rec_lapack10,'--','LineWidth', 2, 'color', 'r'); hold on;
   loglog(threads,rec_lapack50,'--','LineWidth', 2, 'color', 'r'); hold on;
   loglog(threads,rec_lapack100,'--','LineWidth', 2, 'color', 'r'); hold on;
   loglog(threads,rec_lapack200,'--','LineWidth', 2, 'color', 'b'); hold on;
   loglog(threads,rec_lapack300,'--','LineWidth', 2, 'color', 'g'); hold on;
   loglog(threads,rec_lapack400,'--','LineWidth', 2, 'color', 'c'); hold on;
   loglog(threads,rec_lapack500,'--','LineWidth', 2, 'color', 'm'); hold on;
   loglog(threads,rec_lapack600,'-*','LineWidth', 2, 'color', 'r'); hold on;
   loglog(threads,rec_lapack700,'-*','LineWidth', 2, 'color', 'b'); hold on;
   loglog(threads,rec_lapack800,'-*','LineWidth', 2, 'color', 'g'); hold on;
   loglog(threads,rec_lapack900,'-*','LineWidth', 2, 'color', 'c'); hold on;
   loglog(threads,rec_lapack1k,'-*','LineWidth', 2, 'color', 'm'); hold on;
   loglog(threads,rec_lapack2k,'--','LineWidth', 2, 'color', 'y'); hold on;
   loglog(threads,rec_lapack5k,'--','LineWidth', 2, 'color', 'y'); hold on;
   loglog(threads,rec_lapack10k,'--','LineWidth', 2, 'color', 'y'); hold on;
   grid on;
   axis([ 1 16 0 1000 ]);
   legend('MKL LAPACK', 'mt = 1', 'mt = 10', 'mt = 50', 'mt = 100', 'mt = 200', 'mt = 300', 'mt = 400', 'mt = 500', 'mt = 600', 'mt = 700', 'mt = 800', 'mt = 900', 'mt = 1000', 'mt = 2500', 'mt = 5000', 'mt = 10000');
   title({'Profiling of QR Factorization Using LAPACK','Timing Results','Comparing MKL LAPACK to Recursive LAPACK'});
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
   loglog(threads,mkl_lapack,'LineWidth', 2, 'color', 'k'); hold on;
   loglog(threads,rec_lapack1,'--','LineWidth', 2, 'color', 'b'); hold on;
   loglog(threads,rec_lapack10,'--','LineWidth', 2, 'color', 'g'); hold on;
   loglog(threads,rec_lapack100,'--','LineWidth', 2, 'color', 'r'); hold on;
   loglog(threads,rec_lapack500,'--','LineWidth', 2, 'color', 'm'); hold on;
   loglog(threads,rec_lapack1k,'--','LineWidth', 2, 'color', 'c'); hold on;
   loglog(threads,rec_lapack10k,'--','LineWidth', 2, 'color', 'y'); hold on;
   grid on;
   axis([ 1 16 0 1000 ]);
   legend('MKL LAPACK', 'mt = 1', 'mt = 10', 'mt = 100', 'mt = 500', 'mt = 1000', 'mt = 10000');
   title({'Profiling of QR Factorization Using LAPACK','Timing Results','Comparing MKL LAPACK to Recursive LAPACK'});
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
   plot(threads,rec_lapack_gflop1,'--','LineWidth', 2, 'color', 'r'); hold on;
   plot(threads,rec_lapack_gflop10,'--','LineWidth', 2, 'color', 'r'); hold on;
   plot(threads,rec_lapack_gflop50,'--','LineWidth', 2, 'color', 'r'); hold on;
   plot(threads,rec_lapack_gflop100,'--','LineWidth', 2, 'color', 'r'); hold on;
   plot(threads,rec_lapack_gflop200,'--','LineWidth', 2, 'color', 'b'); hold on;
   plot(threads,rec_lapack_gflop300,'--','LineWidth', 2, 'color', 'g'); hold on;
   plot(threads,rec_lapack_gflop400,'--','LineWidth', 2, 'color', 'c'); hold on;
   plot(threads,rec_lapack_gflop500,'--','LineWidth', 2, 'color', 'm'); hold on;
   plot(threads,rec_lapack_gflop600,'-*','LineWidth', 2, 'color', 'r'); hold on;
   plot(threads,rec_lapack_gflop700,'-*','LineWidth', 2, 'color', 'b'); hold on;
   plot(threads,rec_lapack_gflop800,'-*','LineWidth', 2, 'color', 'g'); hold on;
   plot(threads,rec_lapack_gflop900,'-*','LineWidth', 2, 'color', 'c'); hold on;
   plot(threads,rec_lapack_gflop1k,'-*','LineWidth', 2, 'color', 'm'); hold on;
   plot(threads,rec_lapack_gflop2k,'--','LineWidth', 2, 'color', 'y'); hold on;
   plot(threads,rec_lapack_gflop5k,'--','LineWidth', 2, 'color', 'y'); hold on;
   plot(threads,rec_lapack_gflop10k,'--','LineWidth', 2, 'color', 'y'); hold on;
   grid on;
   axis([ 1 16 0 230 ]);
   legend('MKL LAPACK', 'mt = 1', 'mt = 10', 'mt = 50', 'mt = 100', 'mt = 200', 'mt = 300', 'mt = 400', 'mt = 500', 'mt = 600', 'mt = 700', 'mt = 800', 'mt = 900', 'mt = 1000', 'mt = 2500', 'mt = 5000', 'mt = 10000');
   title({'Performance Using One Node','m = n = 10,000, varying mt and number of threads'});
   xlabel('number of threads')
   ylabel('performance for varying mt sizes - (GFlop/sec)')
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
   plot(threads,rec_lapack_gflop1,'--','LineWidth', 2, 'color', 'b'); hold on;
   plot(threads,rec_lapack_gflop10,'--','LineWidth', 2, 'color', 'g'); hold on;
   plot(threads,rec_lapack_gflop100,'--','LineWidth', 2, 'color', 'r'); hold on;
   plot(threads,rec_lapack_gflop500,'--','LineWidth', 2, 'color', 'm'); hold on;
   plot(threads,rec_lapack_gflop1k,'--','LineWidth', 2, 'color', 'c'); hold on;
   plot(threads,rec_lapack_gflop10k,'--','LineWidth', 2, 'color', 'y'); hold on;
   grid on;
   axis([ 1 16 0 230 ]);
   legend('MKL LAPACK', 'mt = 1', 'mt = 10', 'mt = 100', 'mt = 500', 'mt = 1000', 'mt = 10000');
   title({'Performance Using One Node','m = n = 10,000, varying mt and number of threads'});
   xlabel('number of threads')
   ylabel('performance for varying mt sizes - (GFlop/sec)')
   set( ax, 'XTick', threads );
   set( ax, 'XTickLabel', threads );
   set(gca,'FontSize',12);
   grid on

   return
