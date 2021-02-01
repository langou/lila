%
   close all
   clear
%
%  This is a comparison of recursive DGEQRF and MKL DGEQRF
%
  threads = [ 1, 2, 4, 8, 16 ];
  mt = 100;
  m = [ 1000:1000:10000];
  n = m;
%
%  
%
   mkl_lapack1  = [ 0.1653, 1.2117, 3.8489, 8.8892, 16.9054, 28.7717, 44.9727, 66.9425, 93.4293, 128.9239 ];
   rec_lapack1  = [ 0.1625, 1.0875, 3.5581, 8.3382, 16.0435, 27.3657, 43.3848, 64.7145, 91.1293, 124.6622 ];
%   
   mkl_lapack2  = [ 0.1069, 0.6534, 2.0758, 4.6731,  8.7800, 15.0983, 24.1928, 35.4736, 48.6494,  65.7663 ];
   rec_lapack2  = [ 0.1254, 0.7120, 2.2853, 5.1460,  9.5088, 15.8336, 24.4775, 35.5100, 49.9400,  67.4299 ];
%   
   mkl_lapack4  = [ 0.0653, 0.3761, 1.1525, 2.6167,  4.7184,  8.1440, 12.4787, 18.5740, 26.5122,  35.0425 ];
   rec_lapack4  = [ 0.0902, 0.4450, 1.3467, 2.6623,  5.0609,  8.3391, 13.0657, 19.3204, 27.1474,  36.3665 ];
%   
   mkl_lapack8  = [ 0.0445, 0.2396, 0.6653, 1.4776,  2.8689,  4.6477,  7.0075, 10.4828, 14.2883,  19.6928 ];
   rec_lapack8  = [ 0.0712, 0.3053, 0.8233, 1.7416,  3.1400,  5.2313,  7.7375, 11.0440, 15.2604,  20.6116 ];
%   
   mkl_lapack16 = [ 0.0375, 0.1776, 0.4716, 0.9846,  1.8478,  3.0658,  4.5590,  6.4314,  8.8718,  11.7264 ];
   rec_lapack16 = [ 0.0637, 0.2613, 0.6933, 1.2615,  2.2793,  3.7106,  5.2486,  7.2265, 10.1446,  13.4670 ];
%   
   mkl_lapack_gflop1  = [ 16.1252,  17.6066,  18.7069,  19.1994,  19.7176,  20.0197,  20.3383,  20.3956,  20.8072,  20.6840 ];
   rec_lapack_gflop1  = [ 16.4081,  19.6163,  20.2354,  20.4680,  20.7769,  21.0482,  21.0827,  21.0978,  21.3323,  21.3911 ];
%
   mkl_lapack_gflop2  = [ 24.9468,  32.6486,  34.6851,  36.5212,  37.9649,  38.1501,  37.8074,  38.4888,  39.9594,  40.5476 ];
   rec_lapack_gflop2  = [ 21.2673,  29.9647,  31.5062,  33.1638,  35.0552,  36.3783,  37.3677,  38.4493,  38.9267,  39.5472 ];
%
   mkl_lapack_gflop4  = [ 40.0767,  56.7223,  62.4757,  65.2213,  70.6453,  70.7270,  73.2981,  73.5076,  73.3248,  76.0981 ];
   rec_lapack_gflop4  = [ 29.5583,  47.9447,  53.4652,  64.1043,  65.8651,  69.0722,  70.0052,  70.6679,  71.6091,  73.3275 ];
%
   mkl_lapack_gflop8  = [ 59.8768,  89.0533, 108.2223, 115.5054, 116.1904, 123.9324, 130.5263, 130.2448, 136.0555, 135.4135 ];
   rec_lapack_gflop8  = [ 37.4732,  69.8737,  87.4522,  97.9946, 106.1557, 110.1067, 118.2115, 123.6269, 127.3884, 129.3769 ];
%
   mkl_lapack_gflop16 = [ 71.1280, 120.1181, 152.6834, 173.3396, 180.3952, 187.8794, 200.6306, 212.2903, 219.1219, 227.4065 ];
   rec_lapack_gflop16 = [ 41.8412,  81.6449, 103.8540, 135.2846, 146.2429, 155.2330, 174.2698, 188.9351, 191.6297, 198.0147 ];
%
%
%
   figure
   ax = axes;
   plot(m,mkl_lapack1,'LineWidth', 2, 'color', 'k'); hold on;
   plot(m,rec_lapack1,'--','LineWidth', 2, 'color', 'k'); hold on;
   plot(m,mkl_lapack2,'LineWidth', 2, 'color', 'g'); hold on;
   plot(m,rec_lapack2,'--','LineWidth', 2, 'color', 'g'); hold on;
   plot(m,mkl_lapack4,'LineWidth', 2, 'color', 'r'); hold on;
   plot(m,rec_lapack4,'--','LineWidth', 2, 'color', 'r'); hold on;
   plot(m,mkl_lapack8,'LineWidth', 2, 'color', 'c'); hold on;
   plot(m,rec_lapack8,'--','LineWidth', 2, 'color', 'c'); hold on;
   plot(m,mkl_lapack16,'LineWidth', 2, 'color', 'm'); hold on;
   plot(m,rec_lapack16,'--','LineWidth', 2, 'color', 'm'); hold on;
   grid on;
   axis([ 1000 10000 0 130 ]);
   legend('MKL LAPACK 1 thread','Recursive LAPACK 1 thread', 'MKL LAPACK 2 thread','Recursive LAPACK 2 thread', 'MKL LAPACK 4 thread','Recursive LAPACK 4 thread', 'MKL LAPACK 8 thread','Recursive LAPACK 8 thread', 'MKL LAPACK 16 thread','Recursive LAPACK 16 thread');
   title({'Profiling of QR Factorization Using LAPACK','Timing Results (nx=1, mt=100)','Comparing MKL LAPACK to Recursive LAPACK'});
   xlabel('matrix size')
   ylabel('time (seconds)')
   set( ax, 'XTick', m );
   set( ax, 'XTickLabel', m );
   set(gca,'FontSize',12);
   grid on
%
%
%
   figure
   ax = axes;
   plot(m,mkl_lapack_gflop1,'LineWidth', 2, 'color', 'k'); hold on;
   plot(m,rec_lapack_gflop1,'--','LineWidth', 2, 'color', 'k'); hold on;
   plot(m,mkl_lapack_gflop2,'LineWidth', 2, 'color', 'g'); hold on;
   plot(m,rec_lapack_gflop2,'--','LineWidth', 2, 'color', 'g'); hold on;
   plot(m,mkl_lapack_gflop4,'LineWidth', 2, 'color', 'r'); hold on;
   plot(m,rec_lapack_gflop4,'--','LineWidth', 2, 'color', 'r'); hold on;
   plot(m,mkl_lapack_gflop8,'LineWidth', 2, 'color', 'c'); hold on;
   plot(m,rec_lapack_gflop8,'--','LineWidth', 2, 'color', 'c'); hold on;
   plot(m,mkl_lapack_gflop16,'LineWidth', 2, 'color', 'm'); hold on;
   plot(m,rec_lapack_gflop16,'--','LineWidth', 2, 'color', 'm'); hold on;
   grid on;
   axis([ 1000 10000 0 230 ]);
   legend('MKL LAPACK 1 thread','Recursive LAPACK 1 thread', 'MKL LAPACK 2 thread','Recursive LAPACK 2 thread', 'MKL LAPACK 4 thread','Recursive LAPACK 4 thread', 'MKL LAPACK 8 thread','Recursive LAPACK 8 thread', 'MKL LAPACK 16 thread','Recursive LAPACK 16 thread');
   title({'Performance Using One Node','nx = 1, mt = 100; varying matrix size'});
   xlabel('matrix size')
   ylabel('performance - (GFlop/sec)')
   set( ax, 'XTick', m );
   set( ax, 'XTickLabel', m );
   set(gca,'FontSize',12);
   grid on
%
%
%
   return
