% ---------------------------------------------------------------------------- %
% DESCRIPTION                                                                  %
% -----------                                                                  %
% Returns an average curve of multiple curves by linear interpolation. The     %
% returned average curve also has unique and sorted abscissae.                 %
% ---------------------------------------------------------------------------- %
% INPUT                                                                        %
% -----                                                                        %
% abscissae:                                                                   %
% If this argument is a scalar, it specifies the number of abscissae on which  %
% the average curve is evaluated, automatically ranging from the global        %
% minimum to the global maximum abscissa considering all curves. If this       %
% argument is a vector, the average curve is evaluated in the abscissae        %
% specified by it.                                                             %
%                                                                              %
% curves:                                                                      %
% A cell array of length n containing the data of each curve, where n is the   %
% number of curves. The data of each curve is defined by a two-column matrix,  %
% where the first column corresponds to the abscissae and the second column    %
% corresponds to the ordinates. Each curve may have a different number of      %
% points.                                                                      %
% ---------------------------------------------------------------------------- %
% OUTPUT                                                                       %
% ------                                                                       %
% curve:                                                                       %
% The abscissae and ordinates of the average curve in a two-column matrix.     %
% ---------------------------------------------------------------------------- %

function curve = avgcurve(abscissae, curves)

% determine abscissae
if isscalar(abscissae)
    xmin = +Inf;
    xmax = -Inf;
    for i = 1 : length(curves)
        imin = min(curves{i}(:,1));
        imax = max(curves{i}(:,1));
        if imin < xmin, xmin = imin; end
        if imax > xmax, xmax = imax; end
    end
    x = linspace(xmin, xmax, abscissae)';
else
    x = abscissae(:);
end

% compute ordinates
y = zeros(size(x));
k = zeros(size(x));
for i = 1 : length(curves)
    [xi, ii] = unique(curves{i}(:,1), "first");
    yi = curves{i}(ii,2);
    interp = interp1(xi, yi, x, "linear");
    k = k + ~isnan(interp);
    interp(isnan(interp)) = 0;
    y = y + interp;
end
y = y ./ k;

% done
curve = [x y];
