n = 5.1; N = 2001;
t = linspace(0, 20*pi, N);  
R = 1; r = R/n;
x1 = (2*R-r)*cos(t) + r*cos((n-1)*t); y1 = (2*R-r)*sin(t) - r*sin((n-1)*t);
R = R*(2 - 2/n); r = R/n;
x2 = (R-2*r)*cos(t) + 2*r*cos((n-1)*t); y2 = (R-2*r)*sin(t) - 2*r*sin((n-1)*t);
% R = R*(1 - 2/n); 
% n = 3.1; r = R/n;
% x3 = (2*R-r)*cos(t) + r*cos((n-1)*t); y3 = (2*R-r)*sin(t) - r*sin((n-1)*t);
figure(Color='w'); 
plt1 = plot(x1(1), y1(1), 'b', LineWidth = 2);  hold on; axis equal;
plt2 = plot(x2(1), y2(1), 'r', LineWidth = 2);
% plt3 = plot(x3(1), y3(1), 'k', LineWidth = 2);
axis(2*[-1,1,-1,1]); 
for i = 2:N
    plt1.XData = x1(1:i); plt1.YData = y1(1:i);
    plt2.XData = x2(1:i); plt2.YData = y2(1:i);
%     plt3.XData = x3(1:i); plt3.YData = y3(1:i);
    drawnow
end