function Pen = MakePen(PaperLevel, xcoord, ycoord, Hg)
[penx, peny, penz] = cylinder([0,0.3,0.3,0], 8); 
penz(1,:) = 0; penz(2,:) = 1; penz(3:4,:) = 7; 
pen = surf(xcoord + penx, ycoord + peny, PaperLevel + penz+0.1, 'FaceColor', 'r');
set(pen,'Parent', Hg); Pen.pen = pen;

Pen.penup = @(zup) zupp(zup);

t1 = (90:20:270)*2*pi/360;
y = [5, 0.5*cos(t1), 5, 5]; x = [0.5, 0.5*sin(t1), -0.5, 0.5];
X = xcoord + [x;x]; Y = ycoord + [y;y]; Z = 1+[0;0.2]*ones(1,numel(y));
penholder = surf(X, Y, Z); penholder.FaceColor = 0.3*[1,1,1];
PHend1 = patch(X(1, :), Y(1, :), Z(1, :), 0.7*[1,1,1]); 
PHend2 = patch(X(2, :), Y(2, :), Z(2, :), 0.7*[1,1,1]); 
set(penholder,'Parent', Hg); set(PHend1,'Parent', Hg); set(PHend2,'Parent', Hg); 
Pen.penholder = penholder;

    function zupp(zup)
        pen.ZData = pen.ZData + zup;
    end
end

