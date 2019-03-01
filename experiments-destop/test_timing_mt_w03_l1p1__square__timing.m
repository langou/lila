   clear
   close all
   X = load('test_timing_w03_l1p1__square__1thread.dat');
   m = X(:,2);
   perf = X(:,6);
   Y = load('test_timing_w03_l1p1__square__2threads.dat');
   m = Y(:,2);
   perf = Y(:,6);
   T = load('test_timing_w03_l1p1__square__3threads.dat');
   m = T(:,2);
   perf = T(:,6);
   Z = load('test_timing_w03_l1p1__square__4threads.dat');
   m = Z(:,2);
   perf = Z(:,6);

   ax = axes;
   loglog( [1 2 3 4 ], [ X(15,6) Y(15,6) T(15,6) Z(15,6) ], 'r', 'LineWidth', 2); hold on;
   %loglog( [1 2 3 4 ], [ X(15,6) X(15,6)/2 X(15,6)/3 X(15,6)/4 ], 'r--', 'LineWidth', 2); hold on;
   loglog( [1 2 3 4 ], [ X(13,6) Y(13,6) T(13,6) Z(13,6) ], 'b', 'LineWidth', 2); hold on;
   loglog( [1 2 3 4 ], [ X(10,6) Y(10,6) T(10,6) Z(10,6) ], 'g', 'LineWidth', 2); hold on;
   loglog( [1 2 3 4 ], [ X(5,6) Y(5,6) T(5,6) Z(5,6) ], 'm', 'LineWidth', 2); hold on;
   loglog( [1 2 3 4 ], [ X(2,6) Y(2,6) T(2,6) Z(2,6) ], 'k', 'LineWidth', 2); hold on;
   %loglog( [1 2 3 4 ], [ X(1,6) Y(1,6) T(1,6) Z(1,6) ], 'm', 'LineWidth', 2); hold on;
   %loglog( [1 2 3 4 ], [ X(1,6) X(1,6)/2 X(1,6)/3 X(1,6)/4 ], 'm--', 'LineWidth', 2); hold on;
   grid on
   set( ax, 'XTick',[ 1 2 3 4]);
   axis([1 4 0 45])

   legend
   legend('mt = 5000','mt = 3000','mt = 1000','mt = 500','mt = 200')
   title('scalability experiment, m=n=5000, VRQT, mt')
   xlabel('number of processors')
   ylabel('time in seconds')

   figure
   ax = axes;
   loglog( [1 2 3 4 ], [ X(15,7) Y(15,7) T(15,7) Z(15,7) ], 'r', 'LineWidth', 2); hold on;
   %loglog( [1 2 3 4 ], [ X(15,7) 2*X(15,7) 3*X(15,7) 4*X(15,7) ], 'r--', 'LineWidth', 2); hold on;
   loglog( [1 2 3 4 ], [ X(13,7) Y(13,7) T(13,7) Z(13,7) ], 'b', 'LineWidth', 2); hold on;
   loglog( [1 2 3 4 ], [ X(10,7) Y(10,7) T(10,7) Z(10,7) ], 'g', 'LineWidth', 2); hold on;
   loglog( [1 2 3 4 ], [ X(5,7) Y(5,7) T(5,7) Z(5,7) ], 'm', 'LineWidth', 2); hold on;
   loglog( [1 2 3 4 ], [ X(2,7) Y(2,7) T(2,7) Z(2,7) ], 'k', 'LineWidth', 2); hold on;
   %loglog( [1 2 3 4 ], [ X(1,7) Y(1,7) T(1,7) Z(1,7) ], 'r', 'LineWidth', 2); hold on;
   %loglog( [1 2 3 4 ], [ X(1,7) X(1,7)/2 X(1,7)/3 X(1,7)/4 ], 'r--', 'LineWidth', 2); hold on;
   grid on
   set( ax, 'XTick',[ 1 2 3 4]);
   axis([1 4 0 45])

   legend
   legend('mt = 5000','mt = 3000','mt = 1000','mt = 500','mt = 200')
   title('scalability experiment, m=n=5000, VRQT, mt')
   xlabel('number of processors')
   ylabel('GFlop/second')



   figure
   ax = axes;
   loglog( X(1:15,3), X(1:15,7) , 'r', 'LineWidth', 2); hold on;
   loglog( X(1:15,3), Y(1:15,7) , 'g', 'LineWidth', 2); hold on;
   loglog( X(1:15,3), T(1:15,7) , 'c', 'LineWidth', 2); hold on;
   loglog( X(1:15,3), Z(1:15,7) , 'b', 'LineWidth', 2); hold on;
   grid on
   set( ax, 'XTick',[ 1 200 500 1000 3000 5000]);
   axis([1 5000 0 35])

   legend
   legend('1 core','2 cores','3 cores','4 cores')
   title('scalability experiment, m=n=5000, VRQT, mt')
   xlabel('mt')
   ylabel('GFlop/second')

   figure
   ax = axes;
   loglog( X(1:15,3), X(1:15,6) , 'r', 'LineWidth', 2); hold on;
   loglog( X(1:15,3), Y(1:15,6) , 'g', 'LineWidth', 2); hold on;
   loglog( X(1:15,3), T(1:15,6) , 'c', 'LineWidth', 2); hold on;
   loglog( X(1:15,3), Z(1:15,6) , 'b', 'LineWidth', 2); hold on;
   grid on
   set( ax, 'XTick',[ 1 200 500 1000 3000 5000]);
   axis([1 5000 0 350])

   legend
   legend('1 core','2 cores','3 cores','4 cores')
   title('scalability experiment, m=n=5000, VRQT, mt')
   xlabel('mt')
   ylabel('GFlop/second')
