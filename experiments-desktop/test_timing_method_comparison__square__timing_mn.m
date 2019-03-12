   clear
   close all
   X = load('test_timing_w03_l1p1_mn__square__1thread.dat');
   Y = load('test_timing_lapack__square__1thread.dat');
   Z = load('test_timing_w03_qr3__square__1thread.dat');
   ax = axes;
   plot( X(:,1),X(:,6) , 'r', 'LineWidth', 2); hold on;
   plot( X(:,1),Y(:,4) , 'g', 'LineWidth', 2); hold on;
   plot( X(:,1),Z(:,3) , 'b', 'LineWidth', 2);
   grid on
   axis([0 5000 0 50])

   legend
   legend('w03 l1p1','LAPACK','qr3')
   title('VRQT, method comparison, m=n')
   xlabel('dimension size')
   ylabel('time in seconds')

   figure
   ax = axes;
   plot( [0;X(:,1)],[0;X(:,7)] , 'r', 'LineWidth', 2); hold on; 
   plot( X(:,1),Y(:,5) , 'g', 'LineWidth', 2); hold on;
   plot( [0;X(:,1)],[0;Z(:,4)] , 'b', 'LineWidth', 2); 
   grid on
   axis([0 5000 0 11])

   legend
   legend('w03 l1p1','LAPACK','qr3')   
   title('VRQT, method comparison, m=n')
   xlabel('dimension size')
   ylabel('GFlop/second')

