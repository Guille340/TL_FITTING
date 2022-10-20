%
%  pfit = polyfitcon(x,y,p0,ip) 
%
%  DESCRIPTION: finds the coefficients of a polynomial PFIT of degree 
%  length(p0) that fits the dataset (X,Y) in a least-squares sense. An
%  initial set of coefficients P0 must be specified. Any number of 
%  polynomial coefficients can be fixed with the logical vector IP. TRUE
%  TRUE elements in IP indicate what coefficient in P0 are fixed (i.e. 
%  PFIT(IP) == P0(IP)). POLYFITCON doesn't support boxed constraints,
%  only constraints in the form of fixed coefficients.
%
%  INPUT VARIABLES
%  - x: vector of data values for the independent variable
%  - y: data values at x
%  - p0: vector of initial polynomial coefficients, in descending order.
%    Y = P0(1)*X^N + P0(2)*X^(N-1) + ... + P0(N)*X + P0(N+1)
%  - ip: logical vector with TRUE values indicating the elements in P0
%    set as fixed polynomial coefficients (PFIT(IP) == P0(IP)).
%
%  OUTPUT VARIABLES
%  - pfit: vector of polynomial coefficients calculated from a least-
%    square fitting of input data (X,Y).
%
%  FUNCTION DEPENDENCIES
%  - None
%
%  TOOLBOX DEPENDENCIES
%  - MATLAB (Core)
%
%  VERSION 1.0
%  Guillermo Jimenez Arranz
%  email: gjarranz@gmail.com
%  28 Apr 2020

function pfit = polyfitcon(x,y,p0,ip) 
    % Ignore NaN values for the Fitting
    ival = ~isnan(y);
    x = x(ival);
    y = y(ival);
   
    % Calculate the polynomial coefficients
    ip = logical(ip); % convert ip into a logical vector
    objFunc = @(peval) sum((y - polyFunc(peval,x,p0,ip)).^2); % objective function (quadratic error)
    p0_eval = p0(~ip); % initial polynomial coefficients (unfixed)
    pfit_eval = fminsearch(objFunc,p0_eval); % polynomial coefficients (unfixed)
    pfit(ip) = p0(ip); % polynomial coefficients to initial values (fixed)
    pfit(~ip) = pfit_eval; % fitted polynomial coefficients (all)
end

function y = polyFunc(peval,x,p0,ip)
    p(ip) = p0(ip); % evaluation polynomial coefficients (fixed)
    p(~ip) = peval; % evaluation polynomial coefficients (unfixed)
    y = polyval(p,x); % evaluation result for values of x
end