function P = fitBezierCurve(x, y)

    n = length(x);
    
    P0 = [x(1); y(1)];
    P1 = [x(floor(n/2)); y(floor(n/2))];
    P2 = [x(end); y(end)];
    
    % Quadratic Bezier curve definition
    B = @(P0,P1,P2,t) ((1-t).^2 .* P0 + 2*(1-t).*t .* P1 + t.^2 .* P2);
    Bx = @(P0,P1,P2,t) ([1,0] * B(P0,P1,P2,t'))';
    By = @(P0,P1,P2,t) ([0,1] * B(P0,P1,P2,t'))';
           
    c0 = [P0; P1; P2; linspace(0, 1, n)'];
    
    error_x = @(c) Bx(c(1:2), c(3:4), c(5:6), c(7:end)) - x;
    error_y = @(c) By(c(1:2), c(3:4), c(5:6), c(7:end)) - y;
    error_t = @(c) c(7:end) - c0(7:end);
    obj = @(c) [error_x(c); error_y(c); error_t(c)];
    
    cSol = lsqnonlin(obj, c0, [-ones(6,1)*inf; zeros(n,1)], [ones(6,1)*inf; ones(n,1)]);
    
    P0 = cSol(1:2);
    P1 = cSol(3:4);
    P2 = cSol(5:6);
        
    P = [P0'; P1'; P2'];
    
end