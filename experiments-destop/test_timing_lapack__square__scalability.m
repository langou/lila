   clear
   close all
   X = load('test_timing_lapack__square__1thread.dat');
   m = X(:,2);
   perf = X(:,5);
   Y = load('test_timing_lapack__square__2threads.dat');
   m = Y(:,2);
   perf = Y(:,5);
   Z = load('test_timing_lapack__square__4threads.dat');
   m = Z(:,2);
   perf = Z(:,5);

   ax = axes;
   loglog( [1 2 4 ], [ X(13,4) Y(13,4) Z(13,4) ], 'r', 'LineWidth', 2); hold on;
   loglog( [1 2 4 ], [ X(13,4) X(13,4)/2 X(13,4)/4 ], 'r--', 'LineWidth', 2); hold on;
   grid on
   set( ax, 'XTick',[ 1 2 3 4]);
   axis([1 4 9 40])

   title('scalability experiment, m=n=5000, VRQT')
   xlabel('number of processors')
   ylabel('time in seconds')

