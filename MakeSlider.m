function Slider = MakeSlider(xcoord, ycoord, Hg)
XY = load("XYcsxn.txt");
XY = XY(4:13,:);
y = XY(:,1)'; z = XY(:,2)'; 
Y = [y;y]+ycoord; Z = [z;z]; X = xcoord + 0.25*[-1;1]*ones(1,numel(y));
baseslider = surf(X, Y, Z); baseslider.FaceColor = 0.3*[1,1,1];
BSend1 = patch(X(1, :), Y(1, :), Z(1, :), 0.3*[1,1,1]); 
BSend2 = patch(X(2, :), Y(2, :), Z(2, :), 0.3*[1,1,1]); 
set(baseslider,'Parent', Hg); set(BSend1,'Parent', Hg); set(BSend2,'Parent', Hg); 
Slider.baseslider = baseslider;