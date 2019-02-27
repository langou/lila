   clear
   close all
   X = load('test_timing_lapack__square__1thread.dat');
   m = X(:,2);
   perf = X(:,5);
   plot( m, perf, 'r', 'LineWidth', 2); hold on;
   X = load('test_timing_lapack__square__2threads.dat');
   m = X(:,2);
   perf = X(:,5);
   plot( m, perf, 'r', 'LineWidth', 2); hold on;
   grid on
   X = load('test_timing_lapack__square__4threads.dat');
   m = X(:,2);
   perf = X(:,5);
   plot( m, perf, 'r', 'LineWidth', 2); hold on;
   grid on
