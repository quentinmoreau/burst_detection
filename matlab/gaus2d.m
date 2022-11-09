function z=gaus2d(x, y, mx, my, sx, sy)
% GAUS2D  Two-dimensional gaussian function
%   x: x grid
%   y: y grid
%   mx: mean in x dimension
%   my: mean in y dimension
%   sx: standard deviation in x dimension
%   sy: standard deviation in y dimension
% returns: Two-dimensional Gaussian distribution
    z=exp(-((x - mx).^2. / (2. * sx^2.) + (y - my).^2. / (2. * sy^2.)));
end
