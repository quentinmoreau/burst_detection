function o=overlap(a,b)
% Find if two ranges overlap
%   a: first range [low, high]
%   b: second range [low, high]
% returns: True if ranges overlap, false otherwise
    o=(a(1)<=b(1) && b(1)<=a(2)) || (b(1)<=a(1) && a(1)<=b(2));
end