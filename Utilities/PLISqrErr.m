function err = SquareErr(coeffs, x1, y1, x2, y2)
    % Interpolation of 'y2' with scaled 'x2' into the domain 'x1' 
    y2sampledInx1 = interp1(1*x2,y2,x1);
    % Squred error calculation
    err = sum((coeffs(2)*y2sampledInx1-y1).^2);
end
