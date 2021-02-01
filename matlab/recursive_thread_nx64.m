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
   rec_lapack10  = [  ];
   rec_lapack100 = [  ];
   rec_lapack500 = [  ];
   rec_lapack1k  = [  ];
   rec_lapack5k  = [  ];
   rec_lapack10k = [  ];

%   
   mkl_lapack_gflop    = [ 20.6840, 40.5476, 76.0981, 135.4135, 227.3878 ];
   rec_lapack_gflop10  = [  ];
   rec_lapack_gflop100 = [  ];
   rec_lapack_gflop500 = [  ];
   rec_lapack_gflop1k  = [  ];
   rec_lapack_gflop5k  = [  ];
   rec_lapack_gflop10k = [  ];
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
