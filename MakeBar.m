function Bar  = MakeBar(S, L, Hg)
XY = load("XYcsxn.txt");
y = XY(:,1)'; z = XY(:,2)';
Y = [y;y]; X = [S*ones(size(y)); (L+1)*ones(size(y))]; Z = [z;z];
bar = surf(X, Y, Z);  bar.FaceColor = 0.7*[1,1,1];
end1 = patch(S*ones(size(y)), y, z, 0.7*[1,1,1]); 
end2 = patch((L+1)*ones(size(y)), y, z, 0.5*[1,1,1]); 
set(bar,'Parent', Hg); set(end1,'Parent', Hg); set(end2,'Parent', Hg);
Bar.bar = bar; Bar.end1 = end1; Bar.end2 = end2;