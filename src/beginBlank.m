function y1 = beginBlank(chunk1)

    x = chunk1;
    dx = gradient(x,1); ddx = del2(x,1);

    T = 1:numel(x);
    cnd1 = dx(T) > 5;
    cnd2 = abs(ddx(T)-dx(T)) < 30;
    cnd3 = abs(ddx(T)-dx(T)) > 2;
    temp = find(cnd1 & cnd2 & cnd3 == 1);
    goodIDX = temp(1);

    trunc_x = x(goodIDX:end);

    y = nt_detrend(trunc_x',5,[],[],[],[],60);
    z = 0;
    yy = y;
    while ~((yy(1) > -25) && (yy(1) < 25))
        yy(1) = [];
        z = z + 1;
    end

    y1 = [zeros(30,1); zeros(goodIDX-1,1); zeros(z,1); yy];

end