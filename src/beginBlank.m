function [y1,dx,post_detrend,goodIDX,total_blank] = beginBlank(chunk)
    % artifact removal based on zero crossings of gradient of smoothed signal 
    [xData, yData] = prepareCurveData([], chunk(1:300));
    % Set up fittype and options.
    ft = fittype( 'smoothingspline' );
    opts = fitoptions( 'Method', 'SmoothingSpline' );
    opts.SmoothingParam = 0.035476098588056344;

    % Fit model to data.
    splineFit = fit(xData,yData,ft,opts);
    xfit = splineFit(xData)';
%     dx = gradient(xfit,1); 
    
    dx = gradient(chunk,1); 
    % Find zero crossings of fitted curve
    zci = @(v) find(v(:).*circshift(v(:), [-1 0]) <= 0);  
    T = 1:numel(xfit);
    zr = zci(dx(T));
    
    if ~isempty(zr(zr < 200))  % must be a zero crossing for first blanking         
        % zero crossing must be before 200 samples 
        goodIDX = min(zr(zr < 200));
    else
        goodIDX = 2;
    end

    % truncate actual chunk based on conditions above 
    trunc_x = chunk(goodIDX:end);

    y = nt_detrend(trunc_x',5,[],[],[],[],60);
    z = 0;
    yy = y;
    post_detrend = [zeros(30,1); zeros(goodIDX-1,1); y]';
    % additional conditions
    while ~((yy(1) > -25) && (yy(1) < 25) && chunk(1)>-6000)
        yy(1) = [];
        chunk(1) = [];
        z = z + 1;
    end
    y1 = [zeros(30,1); zeros(goodIDX-1,1); zeros(z,1); yy]';
    total_blank = 30 + goodIDX-1 + z;
    goodIDX = goodIDX + 30;
end