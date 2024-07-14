function [width, height] = fwhm(data, srate)

% the width and height of the autocorrelation function taken from:
% adapted from https://github.com/matiasim/Critical_Slowing_Epilepsy  and simplifide to a single case
[r, lag] = autocorr(data, srate);
x = lag;
y = r;


half_height = 0.5;
%intersecting line
X = [x(1) x(end)];
Y = [half_height half_height];

P_width = InterX([x(:)'; y(:)'], [X(:)'; Y(:)']);

% modified to also obtain the mid hight intersection
P_height = InterX([y(:)'; x(:)'], [X(:)'; Y(:)']);

if(isempty(P_width))
    width = nan;
else
    width = P_width(1,1);
end

if(isempty(P_height))
    height = nan;
else
    height = P_height(1,1);
end

end