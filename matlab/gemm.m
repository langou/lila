%
   close all
   clear
%
%  This is a comparison of recursive DGEQRF and MKL DGEQRF
%
  threads = [ 1, 2, 4, 8, 16 ];
  m = [ 10000:1000:20000 ];
  n = m;
%
%  
%
   gemm1  = [ 83.2828, 111.3478, 149.5773, 182.4914, 227.2768, 281.1742, 344.8917, 415.0672, 487.6511, 569.3371, 658.0373 ];
   gemm2  = [ 43.3054,  56.9166,  76.2916,  92.8972, 115.8122, 141.4888, 170.7133, 205.2959, 243.6641, 286.4578, 333.0431 ];
   gemm4  = [ 22.3978,  29.7101,  39.1234,  48.1371,  59.9487,  73.5034,  87.4586, 104.7473, 126.2494, 146.4130, 172.2291 ];
   gemm8  = [ 11.9256,  15.8694,  21.1238,  26.1004,  32.2936,  39.4672,  47.0886,  56.3081,  67.1271,  78.7523,  92.3147 ];
   gemm16 = [  6.9355,   9.1598,  12.0331,  15.0583,  18.4878,  22.6442,  27.5019,  33.1593,  38.2417,  45.3797,  53.8689 ];
%
%
%
   gemm1_gflops  = [  24.0146,  23.9071,  23.1051,  24.0779,  24.1468,  24.0065,  23.7524,  23.6733,  23.9187,  24.0947,  24.3147 ];
   gemm2_gflops  = [  46.1836,  46.7702,  45.2999,  47.2996,  47.3871,  47.7070,  47.9869,  47.8626,  47.8692,  47.8884,  48.0418 ];
   gemm4_gflops  = [  89.2957,  89.5991,  88.3358,  91.2809,  91.5449,  91.8325,  93.6672,  93.8067,  92.3886,  93.6939,  92.8995 ];
   gemm8_gflops  = [ 167.7070, 167.7447, 163.6072, 168.3499, 169.9411, 171.0282, 173.9700, 174.5043, 173.7600, 174.1918, 173.3201 ];
   gemm16_gflops = [ 288.3719, 290.6180, 287.2083, 291.7986, 296.8448, 298.0895, 297.8709, 296.3268, 305.0077, 302.2939, 297.0173 ];
%
%
%
   figure
   ax = axes;
   loglog(m,gemm1,'LineWidth', 2, 'color', 'b'); hold on;
   loglog(m,gemm2,'LineWidth', 2, 'color', 'r'); hold on;
   loglog(m,gemm4,'LineWidth', 2, 'color', 'k'); hold on;
   loglog(m,gemm8,'LineWidth', 2, 'color', 'g'); hold on;
   loglog(m,gemm16,'LineWidth', 2, 'color', 'c'); hold on;
   grid on;
   axis([ 10000 20000 0 700 ]);
   legend('1 thread', '2 thread', '4 thread', '8 thread', '16 thread');
   title({'Profiling of DGEMM','Timing Results','m = n'});
   xlabel('size of matrices A, B, C used in gemm')
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
   plot(m,gemm1_gflops,'LineWidth', 2, 'color', 'b'); hold on;
   plot(m,gemm2_gflops,'LineWidth', 2, 'color', 'r'); hold on;
   plot(m,gemm4_gflops,'LineWidth', 2, 'color', 'k'); hold on;
   plot(m,gemm8_gflops,'LineWidth', 2, 'color', 'g'); hold on;
   plot(m,gemm16_gflops,'LineWidth', 2, 'color', 'c'); hold on;
   grid on;
   axis([ 10000 20000 0 310 ]);
   legend('1 thread', '2 thread', '4 thread', '8 thread', '16 thread');
   title({'Performance Using One Node'});
   xlabel('size of matrices A, B, C used in gemm')
   ylabel('performance for varying matrix sizes (GFlop/sec)')
   set( ax, 'XTick', m );
   set( ax, 'XTickLabel', m );
   set(gca,'FontSize',12);
   grid on

   return
