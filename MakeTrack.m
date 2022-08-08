function Track = MakeTrack(M, Thickness, Hg)
    t1 = (90:270)*2*pi/360; t2 = (-90:90)*2*pi/360; 
    x1 = [cos(t1),  M*2*pi + cos(t2)];  z1 = [sin(t1), sin(t2)];
    x1 = [x1, x1(1)]; z1 = [z1, z1(1)]; 
    L1 = cumsum([0,sqrt(diff(x1).^2 + diff(z1).^2)]);
    x2 = [0.6*cos(t1),  M*2*pi + 0.6*cos(t2)];  z2 = [0.6*sin(t1), 0.6*sin(t2)];
    x2 = [x2, x2(1)]; z2 = [z2, z2(1)]; 
    L2 = cumsum([0,sqrt(diff(x2).^2 + diff(z2).^2)]);
    xm = 0.5*(x1 + x2); zm = 0.5*(z1 + z2);
    Lm = 0.5*(L1 + L2); dl = L1(end)/((2*M+1)*12);
    l1 = 0:dl:L1(end); x1 = interp1(L1, x1, l1); z1 = interp1(L1, z1, l1);
    l2 = interp1(L1, L2, l1); x2 = interp1(L2, x2, l2); z2 = interp1(L2, z2, l2);
    lm = interp1(L1, Lm, l1); xm = interp1(Lm, xm, lm); zm = interp1(Lm, zm, lm);
    N = numel(lm);
    xin = zeros(1,2*N); zin = xin; xout = xin; zout = xin; 
    for n = 1:N
        if(mod(n, 2) == 1)
            xin(2*n-1:2*n) = [xm(n), x2(n)]; zin(2*n-1:2*n) = [zm(n), z2(n)];
        else
            xin(2*n-1:2*n) = [x2(n), xm(n)]; zin(2*n-1:2*n) = [z2(n), zm(n)];
        end
        xout(2*n-1:2*n) = [x1(n), x1(n)]; zout(2*n-1:2*n) = [z1(n), z1(n)];
    end
    yin = -0.5*ones(1, 2*N); yout = 0.5*ones(1, 2*N); 
    
    Xtrack = [xin; xout; xout; xin; xin];
    Ztrack = [zin; zout; zout; zin; zin];
    Ytrack = Thickness*[yin; yin; yout; yout; yin];
    track = surf(Xtrack, Ytrack, Ztrack); 
    track.FaceColor = 0.5*([1,1,1]);
    set(track,'Parent', Hg);
    ltrack = lm;

    Track.track = track;
    Track.Roll = @(theta, hg, trans_axis)RollTrack(theta, hg, trans_axis);

    [Xi, Zi, Yi] = cylinder([0, 0.6, 0.6,0], 20); 
    Yi([1,2],:) = -0.5; Yi([3,4],:) = 0.5; Yi = Thickness*Yi;
    [idler.surf, idler.line] = MakeGear(Xi, Yi, Zi, [], [1,1,1], 0.5);
    [Xip, Zip, Yip] = cylinder(0.1, 8);  Yip = (2.5*Thickness + 0.1)*Yip; 
    Yip = Yip - (2*Thickness);
    Idlerpin = surf(Xip, Yip, Zip, 'FaceColor', [1,1,1]*0.5, 'EdgeAlpha', 0.3);
    set(idler.surf,'Parent', Hg); set(idler.line,'Parent', Hg);
    set(Idlerpin,'Parent', Hg); Track.idler = idler;

    Ns = 12;
    xs = zeros(1, 2*Ns); zs = xs;
    for n = 1:Ns
        if(mod(n, 2) == 1)
            xs(2*n-1:2*n) = [0.8, 0.6]*cos(pi*n/6); 
            zs(2*n-1:2*n) = [0.8, 0.6]*sin(pi*n/6); 
        else
            xs(2*n-1:2*n) = [0.6, 0.8]*cos(pi*n/6); 
            zs(2*n-1:2*n) = [0.6, 0.8]*sin(pi*n/6); 
        end
    end
    xs = [xs, xs(1)]; zs = [zs, zs(1)];
    Xs = M*2*pi + [0*xs;xs;xs;0*xs]; Zs = [0*zs;zs;zs;0*zs]; 
    Ys = [-0.5;-0.5;0.5;0.5]*ones(1, 2*Ns + 1);
    [spoket.surf, spoket.line] = MakeGear(Xs, Ys, Zs, [2,3], [1,1,1], 0.5);
    [Xis, Zis, Yis] = cylinder(0.1, 8);  Yis = (3.5*Thickness + 0.1)*Yis; 
    Yis = Yis - (2*Thickness); Xis =  M*2*pi + Xis;
    Spoketpin = surf(Xis, Yis, Zis, 'FaceColor', [1,1,1]*0.5, 'EdgeAlpha', 0.3);
    set(spoket.surf,'Parent', Hg); set(spoket.line,'Parent', Hg);
    set(Spoketpin,'Parent', Hg); Track.spoket = spoket;


    theta = asin(2/3)*180/pi;
    tmh1 = (90:10:180)*pi/180; tmh2 = linspace(360,360-theta,11)*pi/180;
    tmh3 = linspace(180-theta,360+theta,50)*pi/180;
    tmh4 = linspace(180+theta,180,11)*pi/180;
    tmh5 = (0:10:90)*pi/180;
    xmh = [0.2, -0.2, -0.2+0.2*cos(tmh1), -0.4, -0.8+0.4*cos(tmh2), ...
        0.8*cos(tmh3), 0.8+0.4*cos(tmh4), 0.4, 0.2+0.2*cos(tmh5)];
    zmh = [1.6, 1.6, 1.4+0.2*sin(tmh1), sqrt(0.8), sqrt(0.8)+0.4*sin(tmh2), ...
        0.8*sin(tmh3), sqrt(0.8)+0.4*sin(tmh4), 1.4, 1.4+0.2*sin(tmh5)];
    Xmh = M*2*pi + [xmh; xmh]; Zmh = [zmh; zmh]; Ymh = 2.5+0.2*[1;0]*ones(size(xmh));
    mchnholder1 = surf(Xmh, Ymh, Zmh, 'FaceColor', 0.3*[1,1,1]);
    mchnholder1end1 = patch(Xmh(1, :), Ymh(1, :), Zmh(1, :), 0.7*[1,1,1]); 
    mchnholder1end2 = patch(Xmh(2, :), Ymh(2, :), Zmh(2, :), 0.7*[1,1,1]); 
    set(mchnholder1,'Parent', Hg); set(mchnholder1end1,'Parent', Hg); 
    set(mchnholder1end2,'Parent', Hg); 

    Ymh = 1.5+0.2*[1;0]*ones(size(xmh));
    mchnholder2 = surf(Xmh, Ymh, Zmh, 'FaceColor', 0.3*[1,1,1], 'EdgeColor','none');
    mchnholder2end1 = patch(Xmh(1, :), Ymh(1, :), Zmh(1, :), 0.7*[1,1,1]); 
    mchnholder2end2 = patch(Xmh(2, :), Ymh(2, :), Zmh(2, :), 0.7*[1,1,1]); 
    set(mchnholder2,'Parent', Hg); set(mchnholder2end1,'Parent', Hg);
    set(mchnholder2end2,'Parent', Hg); 

    [Xmch, Zmch, Ymch] = cylinder([0,0.6,0.6,0], 20); 
    Ymch([1,2], :) = 0;  Ymch([3,4], :) = 2; Ymch = 1 + Ymch;
    [machine.surf, machine.line] = MakeGear(M*2*pi + Xmch, Ymch, Zmch, [], [1,1,1], 0.3);
    set(machine.surf,'Parent', Hg);  set(machine.line,'Parent', Hg); 

    recx = 0.4*[-0.5,-0.5,0.5,0.5,-0.5]; recz = 1.2+0.2*[0,1,1,0,0];
    Xh = M*2*pi + [recx; recx]; Zh = [recz; recz];
    Yh = [-3;3]*ones(size(recx));
    holderbar = surf(Xh, Yh, Zh, 'FaceColor', 0.3*[1,1,1]);
    holderbarend1 = patch(Xh(1, :), Yh(1, :), Zh(1, :), 0.7*[1,1,1]); 
    holderbarend2 = patch(Xh(2, :), Yh(2, :), Zh(2, :), 0.7*[1,1,1]); 
    set(holderbar,'Parent', Hg); set(holderbarend1,'Parent', Hg);
    set(holderbarend2,'Parent', Hg); 

    recx = 0.4*[-1,-0.5,0.5,1,-0.5]; recz = 1+0.2*[0,1,1,0,0];
    Xh = M*2*pi + [recx; recx]; Zh = [recz; recz];
    Yh = [-3;-1]*ones(size(recx));
    baseplatebar = surf(Xh, Yh, Zh, 'FaceColor', 0.5*[1,1,1]);
    baseplatebarend1 = patch(Xh(1, :), Yh(1, :), Zh(1, :), 0.3*[1,1,1]); 
    baseplatebarend2 = patch(Xh(2, :), Yh(2, :), Zh(2, :), 0.3*[1,1,1]); 
    set(baseplatebar,'Parent', Hg); set(baseplatebarend1,'Parent', Hg);
    set(baseplatebarend2,'Parent', Hg); 

    %% Roller
    function RollTrack(theta, hg, trans_axis)
        Matrix  = hg.Matrix;
        if(trans_axis == 'y_axis')
            set(hg,'Matrix', Matrix*makehgtform('translate',[0, -theta, 0]));
        end
        if(trans_axis == 'x_axis')
            set(hg,'Matrix', Matrix*makehgtform('translate',[-theta, 0, 0]));
        end
        ltrack = ltrack + theta;
        %%ltrack = lm + theta;
        ltrack = mod(ltrack, lm(end));
        x1moved = interp1(lm, x1, ltrack); z1moved = interp1(lm, z1, ltrack);
        x2moved = interp1(lm, x2, ltrack); z2moved = interp1(lm, z2, ltrack);
        xmmoved = interp1(lm, xm, ltrack); zmmoved = interp1(lm, zm, ltrack);
        for i = 1:N
            if(mod(i, 2) == 1)
                xin(2*i-1:2*i) = [xmmoved(i), x2moved(i)];
                zin(2*i-1:2*i) = [zmmoved(i), z2moved(i)];
            else
                xin(2*i-1:2*i) = [x2moved(i), xmmoved(i)];
                zin(2*i-1:2*i) = [z2moved(i), zmmoved(i)];
            end
            xout(2*i-1:2*i) = [x1moved(i), x1moved(i)];
            zout(2*i-1:2*i) = [z1moved(i), z1moved(i)];
        end
        track.XData = [xin; xout; xout; xin; xin];
        track.ZData = [zin; zout; zout; zin; zin];

        RotateGear(spoket.surf, spoket.line,theta,[2,3]);
    end

    %% 3d Rotation Matrix
    function M = Mxyz(U,theta)
    c = cos(theta); s = sin(theta); ux = U(1); uy = U(2); uz = U(3);
    M = [ux*ux*(1 - c) +    c  uy*ux*(1 - c) - uz*s  uz*ux*(1 - c) + uy*s
         uy*ux*(1 - c) + uz*s  uy*uy*(1 - c) +    c  uz*uy*(1 - c) - ux*s
         uz*ux*(1 - c) - uy*s  uy*uz*(1 - c) + ux*s  uz*uz*(1 - c) +    c];  
    end

    %% RotateGear
    function RotateGear(gear_surf, gear_line, theta, TeethLine)
        X = gear_surf.XData; Y = gear_surf.YData; Z = gear_surf.ZData; 
        [Nl,Lg1] = size(X); xc = X(1,1); yc = Y(1,1); zc = Z(1,1);
        dx = X(1,1)-X(Nl,1); dy = Y(1,1)-Y(Nl,1); dz = Z(1,1)-Z(Nl,1);
        U = [dx, dy, dz]; U = U/norm(U);
        Mat = Mxyz(U,theta); XYZ = Mat*[X(:)'-xc;Y(:)'-yc;Z(:)'-zc];
        X = reshape(XYZ(1,:), Nl,Lg1)+xc; 
        Y = reshape(XYZ(2,:), Nl,Lg1)+yc; 
        Z = reshape(XYZ(3,:), Nl,Lg1)+zc;
        gear_surf.XData = X; gear_surf.YData = Y; gear_surf.ZData = Z; 
    
        Xl = [X(TeethLine,:);nan(1,Lg1)]; Yl = [Y(TeethLine,:);nan(1,Lg1)]; 
        Zl = [Z(TeethLine,:);nan(1,Lg1)]; Xl = Xl(:); Yl = Yl(:); Zl = Zl(:); 
        for i = 2:Nl-1
            Xl = [Xl; nan; X(i,:)']; Yl = [Yl; nan; Y(i,:)']; Zl = [Zl; nan; Z(i,:)'];
        end
        gear_line.XData = Xl; gear_line.YData = Yl; gear_line.ZData = Zl; 
    end

    %% GearPlotter
    function [gear_surf, gear_line] = MakeGear(X, Y, Z, TeethLine, color, alpha)
        gear_surf = surf(X, Y, Z,'EdgeAlpha',0.,'FaceColor', alpha*color); hold on;
        [Nl,Lg1] = size(X); 
        Xl = [X(TeethLine,:);nan(1,Lg1)]; Yl = [Y(TeethLine,:);nan(1,Lg1)]; 
        Zl = [Z(TeethLine,:);nan(1,Lg1)]; Xl = Xl(:); Yl = Yl(:); Zl = Zl(:); 
        for i = 2:Nl-1
            Xl = [Xl; nan; X(i,:)']; Yl = [Yl; nan; Y(i,:)']; Zl = [Zl; nan; Z(i,:)'];
        end
        gear_line = plot3(Xl, Yl, Zl, 'k'); daspect([1,1,1]);
    end
end





