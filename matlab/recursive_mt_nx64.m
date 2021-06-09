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
   rec_lapack1  = [  ];
%
   mkl_lapack2  = 65.7663*ones(1,size(mt,2));
   rec_lapack2  = [  ];
%
   mkl_lapack4  = 35.0425*ones(1,size(mt,2));
   rec_lapack4  = [  ];
%
   mkl_lapack8  = 19.6928*ones(1,size(mt,2));
   rec_lapack8  = [  ];
%
   mkl_lapack16 = 11.7264*ones(1,size(mt,2));
   rec_lapack16 = [  ];
%   
   mkl_lapack_gflop1  = 20.6840*ones(1,size(mt,2));
   rec_lapack_gflop1  = [  ];
%   
   mkl_lapack_gflop2  = 40.5476*ones(1,size(mt,2));
   rec_lapack_gflop2  = [  ];
%   
   mkl_lapack_gflop4  = 76.0981*ones(1,size(mt,2));
   rec_lapack_gflop4  = [  ];
%   
   mkl_lapack_gflop8  = 135.4135*ones(1,size(mt,2));
   rec_lapack_gflop8  = [  ];
%   
   mkl_lapack_gflop16 = 227.3878*ones(1,size(mt,2));
   rec_lapack_gflop16 = [  ];
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
