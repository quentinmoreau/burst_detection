function [right_loc, left_loc, up_loc, down_loc]=fwhm_burst_norm(tf, peak)
% FWHM_BURST_NORM  Find two-dimensional FWHM
%   tf: TF spectrum
%   peak: peak of activity [freq, time]
% returns: right, left, up, down limits for FWM

    half=tf(peak(1),peak(2))/2;
    
    right_loc = NaN;
    % Find right limit (values to right of peak less than half value at
    % peak)
    cand=find(tf(peak(1),peak(2):end)<=half);
    % If any found, take the first one
    if ~isempty(cand)
        right_loc=cand(1);
    end

    up_loc = NaN;
    % Find up limit (values above peak less than half value at peak)
    cand=find(tf(peak(1):end, peak(2)) <= half);
    % If any found, take the first one
    if ~isempty(cand)
        up_loc=cand(1);
    end

    left_loc = NaN;
    % Find left limit (values below peak less than half value at peak)
    cand=find(tf(peak(1),1:peak(2)-1)<=half);
    % If any found, take the last one
    if ~isempty(cand)
        left_loc = peak(2)-cand(end);
    end

    down_loc = NaN;
    % Find down limit (values below peak less than half value at peak)
    cand=find(tf(1:peak(1)-1,peak(2))<=half);
    % If any found, take the last one
    if ~isempty(cand)
        down_loc = peak(1)-cand(end);
    end

    % Set arms equal if only one found
    if isnan(down_loc)
        down_loc = up_loc;
    end
    if isnan(up_loc)
        up_loc = down_loc;
    end
    if isnan(left_loc)
        left_loc = right_loc;
    end
    if isnan(right_loc)
        right_loc = left_loc;
    end

    % Use the minimum arm in each direction (forces Gaussian to be
    % symmetric in each dimension)
    horiz = min([left_loc, right_loc]);
    vert = min([up_loc, down_loc]);
    right_loc = horiz;
    left_loc = horiz;
    up_loc = vert;
    down_loc = vert;
end