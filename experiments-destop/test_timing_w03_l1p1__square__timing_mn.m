   clear
   close all
   X = load('test_timing_w03_l1p1_mn__square__1thread.dat');
   m = X(:,2);
   perf = X(:,5);

   ax = axes;
   plot( X(:,1),X(:,6) , 'r', 'LineWidth', 2); 
   grid on
   %set( ax, 'XTick',[ 100, 200, 500, X(9:15,1) ]);
   axis([0 5000 0 50])

   title('scalability experiment, VRQT, qr3')
   xlabel('dimension size')
   ylabel('time in seconds')

   figure
   ax = axes;
   plot( X(:,1),X(:,7) , 'r', 'LineWidth', 2); 
   grid on
   %set( ax, 'XTick',[ 100, 200, 500, X(9:15,1) ]);
   axis([0 5000 0 11])

   title('scalability experiment, VRQT, qr3')
   xlabel('dimension size')
   ylabel('GFlop/second')

