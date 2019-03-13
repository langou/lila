   clear
   close all
   X = load('timing_turing_gemm_1thread_openblas.dat');

   figure
   l = line([0 7000],[33.6 33.6]); set(l, 'LineWidth', 2); set(l, 'Color', 'k'); hold on;
   N = X(1:10:end,1);
   n = size(N,1);
   for i = 1:n
      perf(i) = max(X(10*(i-1)+1:10*i,5));
   end
   plot( N, perf, 'r', 'LineWidth', 2); 
   axis([0 5000 0 35])
   grid on

   title('performance of DGEMM, m=n=k, OpenBLAS, 1 thread')
   xlabel('m=n=k')
   ylabel('GFlop/sec')
