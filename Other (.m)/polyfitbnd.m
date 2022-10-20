%  [p,fval] = POLYFITBND(x,y,p0,lb,ub,options)  
%
%  DESCRIPTION: finds the coefficients for a polynomial P(X) of degree 
%  LENGTH(P0) that is a best fit (in a least-squares sense) of the dataset 
%  (X,Y). POLYFITBND is a constrained curve-fitting algorithm; unlike
%  POLYFIT, limits may be set to the polynomial coefficients with input 
%  arguments LB and UB. The function uses FMINSEARCHBND, a multi-
%  dimensional constrained non-linear minimisation algorithm to reach
%  the best-fit solution P. Specific settings for the minimisation
%  algorithm may be applied with the OPTIONS structure (see OPTIMSET,
%  FMINSEARCHBND and FMINSEARCH for details). An initial set of 
%  coefficients P0 must be specified. Fixed polynomial coefficients may
%  be introduced by setting their lower and upper bounds of the same value.
%
%  INPUT VARIABLES
%  - x: vector of data points (independent variable).
%  - y: vector of data values at X.
%  - p0: vector of initial polynomial coefficients, in descending order.
%    Z = P0(1)*Q^N + P0(2)*Q^(N-1) + ... + P0(N)*Q + P0(N+1)
%  - lb: lower bound vector or array, same size as P0. Use -INF for those
%    variables (i.e. elements in X0) with no lower bounds. Leave LB empty 
%    if there are no lower bounds at all. Variables may be fixed by setting
%    the corresponding lower and upper bounds to exactly the same value.
%  - ub: upper bound vector or array, same size as P0. Use INF for those
%    variables (i.e. elements in X0) with no upper bounds. Leave UB empty 
%    if there are no upper bounds at all. Variables may be fixed by setting
%    the corresponding lower and upper bounds to exactly the same value.
%
%  OUTPUT VARIABLES
%  - p: vector of polynomial coefficients that best fits input data (X,Y).
%
%  FUNCTION CALLS
%  p = POLYFITBND(x,y,p0)  
%  p = POLYFITBND(x,y,p0,lb)  
%  p = POLYFITBND(x,y,p0,[],ub)
%  p = POLYFITBND(x,y,p0,lb)  
%  p = POLYFITBND(x,y,p0,lb,ub)  
%  p = POLYFITBND(x,y,p0,lb,ub,options)  
%  [p,fval] = POLYFITBND(x,y,p0,...)  
%
%  FUNCTION DEPENDENCIES
%  - FMINSEARCHBND
%
%  TOOLBOX DEPENDENCIES
%  - MATLAB (Core)
%
%  See also POLYFIT, FMINSEARCHBND

%  VERSION 1.0
%  Guillermo Jimenez Arranz
%  email: gjarranz@gmail.com
%  30 Apr 2020

function [p,fval] = polyfitbnd(x,y,p0,lb,ub,options) 
    % Input Arguments
    if (nargin < 6)
      options = [];
    end
    if nargin < 5
      ub = [];
    end
    if nargin < 4
      lb = [];
    end

    % Ignore NaN values for the Fitting
    ival = ~isnan(y) & ~isnan(x);
    x = x(ival);
    y = y(ival);
   
    % Calculate the polynomial coefficients
    objFunc = @(peval)sum((y - polyval(peval,x)).^2); ...
        % objective function (quadratic error)
    [p,fval] = fminsearchbnd(objFunc,p0,lb,ub,options); % polynomial coefficients (unfixed)
    pause(1)
end