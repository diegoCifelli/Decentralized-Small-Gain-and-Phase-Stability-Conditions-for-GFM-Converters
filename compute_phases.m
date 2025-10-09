function phases = compute_phases(A)
    [rerange, imrange] = nrange_no_plot(A, 500);
    x = rerange + 1i*imrange;
    a = angle(x)';
    angles = mod(a, 2*pi);
    sortedAngles = sort(angles);
    dAngles = diff([sortedAngles; sortedAngles(1) + 2*pi]);
    [maxGap, idx] = max(dAngles);
    minConeAngle = 2*pi - maxGap;
    if idx < length(sortedAngles)
        lowerBound = mod(sortedAngles(idx+1), 2*pi);
    else
        lowerBound = sortedAngles(1);
    end
    upperBound = mod(lowerBound + minConeAngle, 2*pi);
    lowerBound = mod(lowerBound + pi, 2*pi) - pi;
    upperBound = mod(upperBound + pi, 2*pi) - pi;
    phases = [lowerBound upperBound];
    
end