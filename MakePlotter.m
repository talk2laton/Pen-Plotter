function Plotter = MakePlotter
close all;
figure('Position', [1921 41 1920 963], 'Color', 'w');
ax = axes('XLim',[-5 35],'YLim',[-35 30],'ZLim',[-1.2 7]);
view(3); daspect([1,1,1]); hold on; grid on

hg1 = hgtransform('Parent',ax); Track1 = MakeTrack(5, 1, hg1);
hbar1 = hgtransform('Parent',hg1); Bar1 = MakeBar(-3, 5*2*pi, hbar1);
set(hbar1,'Matrix', makehgtform('translate',[0, -2, 0]));
hpen = hgtransform('Parent',hg1); Pen = MakePen(-4, 5*pi/2, -4, hpen);
Slider1 = MakeSlider(5*pi/2, -2,  hpen);
set(hg1,'Matrix', makehgtform('translate',[0, 0, 3]));

hg2 = hgtransform('Parent',ax); Track2 = MakeTrack(9, 1, hg2);
hbar2 = hgtransform('Parent',hg2); Bar2 = MakeBar(-1, 9*2*pi, hbar2);
hsld = hgtransform('Parent',hg1);
Slider2 = MakeSlider(0, 0,  hsld);
set(hsld,'Matrix', makehgtform('translate',[0, -2, -2])* ...
                 makehgtform('zrotate',pi/2));
set(hbar2,'Matrix', makehgtform('translate',[0, -2, 0]));

hbar3 = hgtransform('Parent',hg2); Bar3 = MakeBar(-1, 9*2*pi, hbar3);
set(hbar3,'Matrix', makehgtform('translate',[0, -33, 0]));

HgTrans = makehgtform('zrotate',pi/2)*makehgtform('translate',[-30,2,1]);
set(hg2,'Matrix', HgTrans)


[X, Y] = meshgrid(linspace(2, 29, 28), linspace(-30, 22, 53));
Z = -1+0*X;
surf(X, Y, Z, 'FaceColor', 'w', 'EdgeAlpha',0.3); 

x_current = 7.9;
y_current = -4;
XMin = 4; XMax = 27;
YMin = -28; YMax = 20;
XFac = (XMax - XMin)/(YMax-YMin);

Plotter.plot = @(filename, step, orientation, clr, linewitdth) ...
                 Plot(filename, step, orientation, clr, linewitdth);

    function Plot(filename, step, orientation, clr, linewidth)
        zmove = 0.5; Pen.pen.FaceColor = clr;
        if(ischar(filename))
            [Paths, xmin, xmax, ymin, ymax] = ExtractImagePaths(filename, step);
        else
            [Paths, xmin, xmax, ymin, ymax] = ExtractPaths(filename, step);
        end
        xfac = (xmax-xmin)/(ymax - ymin);
        if(XFac > xfac)
            XMean = 0.5*(XMin + XMax);
            XDel = 0.5*xfac*(YMax - YMin);
            XMin = XMean - XDel;
            XMax = XMean + XDel;
        else
            YMean = 0.5*(YMin + YMax);
            YDel = 0.5*(XMax - XMin)/xfac;
            YMin = YMean - YDel;
            YMax = YMean + YDel;
        end

        if(orientation == "Landscape")
            temp = YMin; YMin = YMax; YMax = temp;
            temp = ymin; ymin = xmin; xmin = temp;
            temp = ymax; ymax = xmax; xmax = temp;
        end
        for n = 1:numel(Paths)
            path = Paths{n}; 
            if(size(path, 1) > 1)
                if(orientation == "Landscape")
                    path = fliplr(path);
                end
                [px, py] = Scale(path(1,1), path(1,2), xmin, ...
                    xmax, ymin, ymax);
                penup(zmove); move(px, py, []); pendn(zmove);
                line3 = plot3(px, py, -1, clr, 'LineWidth', linewidth); 
                for i = 2:size(path,1)
                    [px, py] = Scale(path(i,1), path(i,2), xmin, ...
                        xmax, ymin, ymax); 
                    move(px, py, line3);
                end
            end
        end
        move(4, -28, [])
        if(orientation == "Landscape")
            AZ = linspace(37.5,90);
            EL = linspace(30,90);
            for i = 1:100
                view(AZ(i), EL(i));
                drawnow;
            end
        else
            AZ = linspace(37.5,180);
            EL = linspace(30,90);
            for i = 1:100
                view(AZ(i), EL(i));
                drawnow;
            end
        end
    end
    
    function [px, py] = Scale(xp, yp, xmin, xmax, ymin, ymax)
        fx = (xp-xmin)/(xmax - xmin); fy = (yp-ymin)/(ymax - ymin);
        px = XMin + fx*(XMax - XMin); py = YMin + fy*(YMax - YMin);
    end

    function move(x, y, line3)
        del = sqrt((x_current-x)^2 + (y_current-y)^2);
        n   = ceil(del/0.5);
        if(n > 0)
            dx = (x - x_current)/n; dy = (y - y_current)/n;
            for i = 1:n
                 % Move along x
                 Track1.Roll(-dx, hpen, 'x_axis'); 
                 % move along y
                 Track2.Roll(-dy, hg1, 'y_axis');
                 if(~isempty(line3))
                     line3.XData = [line3.XData, line3.XData(end) + dx];
                     line3.YData = [line3.YData, line3.YData(end) + dy];
                     line3.ZData = [line3.ZData, line3.ZData(end)];
                 end
                 drawnow;
            end
        end
        x_current = x; y_current = y;
    end

    function penup(zmove)
        dz = zmove/10;
        for n = 1:10
            Pen.penup(dz); 
            drawnow;
        end
    end

    function pendn(zmove)
        dz = -zmove/20;
        for n = 1:20
            Pen.penup(dz); 
            drawnow;
        end
    end
end
