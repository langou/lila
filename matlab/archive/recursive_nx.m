%
   close all
   clear
%
%  This is a comparison of recursive DGEQRF and MKL DGEQRF
%
  mt = [ 1, 10, 100, 1000, 5000, 10000 ];
  m = 10000;
  n = 10000;
%
%  
%
   mkl_lapack         = [ 128.0902, 128.0902, 128.0902, 128.0902, 128.0902, 128.0902 ];
   rec_lapack_nx1     = [ 975.0465, 220.1590, 123.8504, 119.1179, 135.2279, 152.0301 ];
   rec_lapack_nx10    = [ , , , , ,  ];
   rec_lapack_nx100   = [ , , 123.4712, 118.6906, 135.2974, 151.8451 ];
   rec_lapack_nx1000  = [ , , 123.2115, 117.6817, 133.7873, 150.4267 ];
   rec_lapack_nx5000  = [ , , 126.0390, 117.5637, 137.2343, 154.1036 ];
   rec_lapack_nx10000 = [ , , 128.5593, 118.7886, 137.8547, 181.4549 ];
%   
   mkl_lapack_gflop         = [ 20.8187, 20.8187, 20.8187, 20.8187, 20.8187, 20.8187 ];
   rec_lapack_gflop_nx1     = [  2.7349, 12.0934, 21.5314, 22.3868, 19.7198, 17.5404 ];
   rec_lapack_gflop_nx10    = [ , , , , ,  ];
   rec_lapack_gflop_nx100   = [ , , 21.5979, 22.4674, 19.7097, 17.5618 ];
   rec_lapack_gflop_nx1000  = [ , , 21.6430, 22.6600, 19.9321, 17.7274 ];
   rec_lapack_gflop_nx5000  = [ , , 21.1575, 22.6827, 19.4315, 17.3044 ];
   rec_lapack_gflop_nx10000 = [ , , 20.7247, 22.4489, 19.3441, 14.6960 ];
%
%
%
   figure
   ax = axes;
   loglog(mt,mkl_lapack,'LineWidth', 2, 'color', 'b'); hold on;
   loglog(mt,rec_lapack_nx1,'LineWidth', 2, 'color', 'r'); hold on;
   grid on;
   axis([ 100 10000 0 153 ]);
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
   plot(mt,mkl_lapack_gflop,'LineWidth', 2, 'color', 'b'); hold on;
   plot(mt,rec_lapack_gflop_nx1,'LineWidth', 2, 'color', 'r'); hold on;
   grid on;
   axis([ 100 10000 0 25 ]);
   legend('MKL LAPACK', 'Recursive - LAPACK Panel');
   title({'Performance Using One Core','m = n = 10,000'});
   xlabel('mt - size of bkock Householder update')
   ylabel('performance for varying mt sizes (GFlop/sec)')
   set( ax, 'XTick', mt );
   set( ax, 'XTickLabel', mt  );
   set(gca,'FontSize',12);
   grid on

   return
