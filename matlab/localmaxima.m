function index = localmaxima(x)
index = find( diff( sign( diff([0; x(:); 0]) ) ) < 0 );