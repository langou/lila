   clear
   close all
   X = load('test_timing_lapack__square__1thread.dat');
   m = X(:,2);
   perf = X(:,5);
   %plot( m, perf, 'r', 'LineWidth', 2); hold on;
   %X = load('test_timing_lapack__square__2threads.dat');
   %m = X(:,2);
   %perf = X(:,5);
   %plot( m, perf, 'r', 'LineWidth', 2); hold on;
   %grid on
   %X = load('test_timing_lapack__square__4threads.dat');
   %m = X(:,2);
   %perf = X(:,5);
   %plot( m, perf, 'r', 'LineWidth', 2); hold on;
   %grid on


   ax = axes;
   plot( X(:,1),X(:,4) , 'r', 'LineWidth', 2); 
   grid on
   %set( ax, 'XTick',[ 100, 200, 500, X(9:15,1) ]);
   axis([0 5000 0 50])

   title('scalability experiment, VRQT, qr3')
   xlabel('dimension size')
   ylabel('time in seconds')

   figure
   ax = axes;
   plot( X(:,1),X(:,5) , 'r', 'LineWidth', 2); 
   grid on
   %set( ax, 'XTick',[ 100, 200, 500, X(9:15,1) ]);
   axis([0 5000 0 10])

   title('scalability experiment, VRQT, qr3')
   xlabel('dimension size')
   ylabel('GFlop/second')
