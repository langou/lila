%
   close all
   clear
%
%  This is a comparison of recursive DGEQRF and MKL DGEQRF
%
  mt = 100;
  m = [ 1000:1000:10000]
  n = m;
%
%  
%
   mkl_lapack1  = [ 0.1653, 1.2117, 3.8489, 8.8892, 16.9054, 28.7717, 44.9727, 66.9425, 93.4293, 128.9239 ];
   rec_lapack1  = [ 0.1625, 1.0875, 3.5581, 8.3382, 16.0435, 27.3657, 43.3848, 64.7145, 91.1293, 124.6622 ];
%   
   mkl_lapack_gflop1  = [ 16.1252,  17.6066,  18.7069,  19.1994,  19.7176,  20.0197,  20.3383,  20.3956,  20.8072,  20.6840 ];
   rec_lapack_gflop1  = [ 16.4081,  19.6163,  20.2354,  20.4680,  20.7769,  21.0482,  21.0827,  21.0978,  21.3323,  21.3911 ];
%
%
%

   figure
   ax = axes;
   plot(m,mkl_lapack1,'LineWidth', 2, 'color', 'k'); hold on;
   plot(m,rec_lapack1,'--','LineWidth', 2, 'color', 'k'); hold on;
   grid on;
   axis([ 1000 10000 0 130 ]);
   legend('MKL LAPACK','Recursive LAPACK');
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
   grid on;
   axis([ 1000 10000 16 22 ]);
   legend('MKL LAPACK', 'Recursive LAPACK');
   title({'Performance Using One Node','nx = 1, mt = 100; varying matrix size'});
   xlabel('matrix size')
   ylabel('performance - (GFlop/sec)')
   set( ax, 'XTick', m );
   set( ax, 'XTickLabel', m );
   set(gca,'FontSize',12);
   grid on

   return
