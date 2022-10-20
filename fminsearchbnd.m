%  [x,fval,exitflag,output] = FMINSEARCHBND(fun,x0,lb,ub,options)
%
%  DESCRIPTION: Multidimensional constrained nonlinear minimization 
%  (Nelder-Mead). Identical to FMINSEARCH, but including bound constraints
%  through transformation. The original constrained fitting variables are 
%  converted into their unconstrained surrogates to be used by FMINSEARCH.
%  Simple continuous functions are used to convert a set of real variables
%  (constrained) to the infinite space of real numbers (unconstrained). 
%  For single and dual bounds, quadratic and sin(X) functions are used, 
%  respectively. FMINSEARCHBND admits fixed variables by simply setting 
%  both bounds (LB and UB) to the exact same value.
%
%  INPUT ARGUMENTS
%  - fun: handle to objective function for which FMINSEARCHBND attempts to 
%    find a local minimiser. FUN accepts an input X and returns a scalar 
%    value FVAL evaluated at X, where X can be a scalar, vector or matrix.
%  - x0: initial set of variables to be evaluated by FUN. X0 can be a
%    scalar, vector or matrix.
%  - lb: lower bound vector or array, same size as X0. Use -INF for those
%    variables (i.e. elements in X0) with no lower bounds. Leave LB empty 
%    if there are no lower bounds at all. Variables may be fixed by setting
%    the corresponding lower and upper bounds to exactly the same value.
%  - ub: upper bound vector or array, same size as X0. Use INF for those
%    variables (i.e. elements in X0) with no upper bounds. Leave UB empty 
%    if there are no upper bounds at all. Variables may be fixed by setting
%    the corresponding lower and upper bounds to exactly the same value.
%  - options: structure with information about the optimisation process. 
%    Created with function OPTIMSET (see OPTIMSET for details). FMINSEARCHBND 
%    uses these options: Display, TolX, TolFun, MaxFunEvals, MaxIter, 
%    FunValCheck, PlotFcns, and OutputFcn.
%
%  OUTPUT ARGUMENTS
%  - x: local minimiser. In other words, a combination of variable values 
%    within the specified bounds (LB,UB) for which FMINSERACHBND found a 
%    local minimum of the objective function FUN.
%  - fval: value of the objective function at the local minimiser X.
%  - exitflag: exit condition. Three possible values of EXITFLAG:
%    ¬  1: maximum coordinate difference between current best point X and 
%      other points in simplex is <= TolX, and corresponding difference in
%      function values is <= TolFun. FMINSEARCHBND stops when it satisfies
%      both TolFun and TolX.
%    ¬  0: maximum number of function evaluations or iterations reached.
%    ¬ -1: algorithm terminated by the output function.
%  - output: structure with the following fields:
%    ¬ iterations: number of iterations taken
%    ¬ funcCount: number of function evaluations
%    ¬ algorithm: algorithm name ('Nelder-Mead simplex direct search')
%    ¬ message: exit message
%
%  FUNCTION CALLS
%   x = FMINSEARCHBND(fun,x0)
%   x = FMINSEARCHBND(fun,x0,lb)
%   x = FMINSEARCHBND(fun,x0,[],ub)
%   x = FMINSEARCHBND(fun,x0,lb,ub)
%   x = FMINSEARCHBND(fun,x0,lb,ub,options)
%  [x,fval] = FMINSEARCHBND(fun,x0,...)
%  [x,fval,exitflag] = FMINSEARCHBND(fun,x0,...)
%  [x,fval,exitflag,output] = FMINSEARCHBND(fun,x0,...)
%
%  FUNCTION DEPENDENCIES
%  - None
%
%  TOOLBOX DEPENDENCIES
%  - MATLAB (Core)
%
%  CONSIDERATIONS & LIMITATIONS
%  - If OPTIONS is supplied, then TolX will apply to the transformed
%    variables. All other FMINSEARCH parameters should be unaffected.
%
%  - Variables which are constrained by both a lower and an upper
%    bound will use a sin transformation. Those constrained by
%    only a lower or an upper bound will use a quadratic
%    transformation, and unconstrained variables will be left alone.
%
%  - Variables may be fixed by setting their respective bounds equal.
%    Any fixed variable will be discarded for solving the problem with
%    FMINSEARCH.
%   
%  - The bounds (LB, UB) are inclusive inequalities, which admit the
%    boundary values themselves, but will not permit ANY function
%    evaluations outside the bounds. These constraints are strictly
%    followed. If your problem has an EXCLUSIVE (strict) constraint which
%    will not admit evaluation at the bound itself, then you must provide
%    a slightly offset bound. An example of this is a function which
%    contains the log of one of its parameters. If you constrain the
%    variable to have a lower bound of zero, then FMINSEARCHBND may
%    try to evaluate the function exactly at zero.
%
%  - Unlike FMINSEARCH, the PROBLEM structure has not been included in
%    FMINSEARCHBND as a call method for simplicity, as it doesn't add much
%    value to the function.
%
%  EXAMPLES
%  The following examples use the Rosenbrock function, a well-known non-
%  convex function used as a performance test problem for optimization 
%  algorithms (Rosenbrock, 1960).
%
%  rosen = @(x) (1-x(1)).^2 + 105*(x(2)-x(1).^2).^2; % Rosenbrock function
%
%  xsol = fminsearchbnd(rosen,[3 3]); % no bounds
%  xsol = fminsearchbnd(rosen,[3 3],[2 2]); lower bound
%  xsol = fminsearchbnd(rosen,[-5 -5],[],[0 0]); % upper bound
%  xsol = fminsearchbnd(rosen,[3 3],[-1 -1],[4 4]); % dual bounds(*)
%  xsol = fminsearchbnd(rosen,[2.5 2.5],[2 2],[3 3]); % dual bounds(**)
%  xsol = fminsearchbnd(rosen,[0 0],[2 2],[3 3]); % dual bounds (***)
%  xsol = fminsearchbnd(rosen,[0 0],[2 -inf],[inf 3]); % mixed bounds
%  xsol = fminsearchbnd(rosen,[3 3],[-1 -1],[-1 4]); % fixed variable
%  [xsol,fval,exitflag,output] = fminsearchbnd(rosen,[3 3],[-1 -1],...
%     [4 4],optimset('Display','iter','TolFun',1.e-12,'TolX',1.e-12,...
%     'FunValCheck','off'));
% 
%  (*) including function minimum
%  (**) not including function minimum
%  (***) with infeasible stgarting guess  
%
%  REFERENCES
%  - Rosenbrock, H.H. (1960). "An automatic method for finding the
%    greatest or least value of a function".
%  - https://uk.mathworks.com/matlabcentral/fileexchange/8277-fminsearchbnd-fminsearchcon
%  - https://stackoverflow.com/questions/56233365/should-i-transform-constraint-optimization-to-unconstrained-optimization
%  - https://math.stackexchange.com/questions/75077/mapping-the-real-line-to-the-unit-interval
%  
%  See also FMINSEARCH, OPTIMSET, FMINBND

%  VERSION 4.1 (20 Apr 2020)
%  - Updated help and code for clarity
%  Guillermo Jimenez Arranz
%  email: gjarranz@gmail.com
%  
%  VERSION 4.0 (23 Jul 2006)
%  John D'Errico
%  email: woodchips@rochester.rr.com

function [x,fval,exitflag,output] = fminsearchbnd(fun,x0,lb,ub,options)

% Input Arguments
xsize = size(x0);
n = numel(x0);
if (nargin < 5) || isempty(options)
  options = optimset('fminsearch');
end
if nargin < 4 || isempty(ub)
  ub = inf(n,1);  
end
if nargin < 3 || isempty(lb)
  lb = -inf(n,1);
end

% Convert Input Variables to Column Vectors
x0 = x0(:);
lb = lb(:);
ub = ub(:);

% Error Control
if ~isa(fun,'function_handle')
    error('Input argument FUN must be a valid function handle')
end
if ~isequal(n,length(lb),length(ub))
  error('X0, LB and UB must all have the same number of elements')
end

% Working parameters structure
params.lb = lb; % lower bound
params.ub = ub; % upper bound
params.fun = fun; % objective function
params.n = n; % number of parameters (varies if one or more are fixed)
params.xsize = xsize;
params.OutputFcn = [];
params.boundType = zeros(n,1);
for i = 1:n
  boundTypeTemp = isfinite(lb(i)) + 2*isfinite(ub(i));
  params.boundType(i) = boundTypeTemp;
  if boundTypeTemp == 3 && lb(i) == ub(i)
    params.boundType(i) = 4;
  end
end
% boundType = 0 for unconstrained variable
% boundType = 1 for lower bound only
% boundType = 2 for upper bound only
% boundType = 3 for dual finite bounds
% boundType = 4 for fixed variable
ieval = params.boundType ~= 4; % index of variables to fit
params.ieval = ieval;

% Transform starting values into their unconstrained surrogates
x0u = directTransform(x0,params); % initial variables (unconstrained)
x0u_eval = x0u(ieval); % initial variables to fit (unconstrained)

if any(ieval) % if one or more variables are not fixed
    % Replace OutputFcn if defined with wrapper output function
    if ~isempty(options.OutputFcn)
      params.OutputFcn = options.OutputFcn;
      outfunAnonymous = @(x)outfun(x,params);
      options.OutputFcn = outfunAnonymous;
    end

    % Call fminsearch, but with our own wrapper input function
    infunAnonymous = @(x)infun(x,params);
    [xu_eval,fval,exitflag,output] = fminsearch(infunAnonymous,...
        x0u_eval,options);
    x = inverseTransform(xu_eval,params); % recover variables in original space
    x = reshape(x,params.xsize);
    
else % if all variables are fixed
    x = inverseTransform(x0u_eval,params); % transform variables to original space
    x = reshape(x,xsize);
    fval = feval(params.fun,x);
    exitflag = 0; % fminsearchbnd was not called
    output.iterations = 0;
    output.funcCount = 1;
    output.algorithm = 'fminsearch';
    output.message = 'All variables were held fixed by the applied bounds';
end

end

% SUBFUNCTIONS
% =====================================================================
function stop = outfun(xu_eval,params)
% DESCRIPTION
% Wrapper for output function. OUTFUN brings the transformed variables 
% XU_EVAL back to the original numerical space, to then evaluate them with
% the user's output function OUTPUTFCN (see OPTIMSET). XU_EVAL includes
% only the variables that need to be evaluated, i.e. those that are
% unfixed (LB ~= UB).
%
% OUTFUN is a nested function and special care must be taken not to use 
% unintentionally variables with identical name to any variable in the 
% parent function FMINSEARCHBND (all variables in the parent function can 
% be accessed in a nested function without the need to pass them as 
% arguments).
%
% INPUT ARGUMENTS
% - xu_eval: vector of polynomial variables in original (constrained) 
%   domain. Only the variables to be solved with FMINSEARCH (i.e. unfixed
%   variables) must be included in XU_EVAL.
%
% OUTPUT ARGUMENTS
% - stop: logical value indicating whether the optimization routine should
%   quit or continue

x = inverseTransform(xu_eval,params); % variables to original domain
stop = params.OutputFcn(x); % call the user supplied OutputFcn  
end

function fval = infun(xu_eval,params)
% DESCRIPTION
% Wrapper for input objective function. INFUN brings the 
% transformed variables XU back to the original numerical space, to
% then evaluated them with the objective function FUN
%
% INPUT ARGUMENTS
% - xu_eval: vector of polynomial variables in original (constrained) 
%   domain. Only the variables to be solved with FMINSEARCH (i.e. unfixed
%   variables) must be included in XU_EVAL.
% - params: structure of working parameters
%
% OUTPUT ARGUMENTS
% - fval: value of the objective function FUN for the given variable
%   combination X (XU_EVAL (unconstrained) -> X (constrained)). 

x = inverseTransform(xu_eval,params); % recover variables in original space
fval = feval(params.fun,reshape(x,params.xsize)); % evaluate FUN

end

function xu = directTransform(x,params)
% DESCRIPTION
% Projects original variables into a transformed numerical space
%
% INPUT ARGUMENTS
% - x: vector of polynomial variables in the original (constrained) domain
% - params: structure of working parameters
%
% OUTPUT ARGUMENTS
% - xu: vector of polynomial variables in transformed (unconstr.) domain

% Load working parameters
n = params.n;
lb = params.lb;
ub = params.ub;
boundType = params.boundType;

% Variables from Original to Transformed Domain
xu = nan(n,1); 
for i = 1:n
  switch boundType(i)
    case 0 % unconstrained variable
      xu(i) = x(i);

    case 1 % lower bound only
      xu(i) = 0; % set to 0 in case that X(i) <= LB(i)
      if x(i) > lb(i)
        xu(i) = sqrt(x(i) - lb(i));
      end

    case 2 % upper bound only
      xu(i) = 0; % set to 0 in case that X(i) >= UB(i)
      if x(i) < ub(i)
        xu(i) = sqrt(ub(i) - x(i));
      end

    case 3 % lower and upper bounds
      if x(i) <= lb(i)
        xu(i) = -pi/2; 
      elseif x(i) >= ub(i)
        xu(i) = pi/2;
      else
        xu(i) = 2*(x(i) - lb(i))/(ub(i)-lb(i)) - 1;
        xu(i) = 2*pi + asin(max(-1,min(1,xu(i)))); 
        % shift by 2*pi to avoid problems at zero in fminsearch. 
        % Otherwise, the initial simplex is vanishingly small
      end
  end
end

end

function x = inverseTransform(xu_eval,params)
% DESCRIPTION
% Brings transformed variables back into the original space
%
% INPUT ARGUMENTS
% - xu_eval: vector of polynomial variables in original (constrained) 
%   domain. Only the variables to be solved with FMINSEARCH (i.e. unfixed
%   variables) must be included in XU_EVAL.
% - params: structure of working parameters
%
% OUTPUT ARGUMENTS
% - x: vector of polynomial variables in original (constrained) domain

% Load working parameters
n = params.n;
xsize = params.xsize;
lb = params.lb;
ub = params.ub;
boundType = params.boundType;
ieval = params.ieval;

% Variables from Transformed to Original Domain
x = nan(xsize);
xu = nan(xsize);
xu(ieval) = xu_eval;
for i = 1:n
  switch boundType(i)
    case 0 % unconstrained variable
        x(i) = xu(i);

    case 1 % lower bound only
        x(i) = lb(i) + xu(i).^2;

    case 2 % upper bound only
        x(i) = params.ub(i) - xu(i).^2;

    case 3 % lower and upper bounds
        x(i) = (sin(xu(i))+1)/2;
        x(i) = x(i)*(ub(i) - lb(i)) + lb(i);
        x(i) = max(lb(i),min(ub(i),x(i))); % avoid float point problems

    case 4 % fixed variable
        x(i) = lb(i); % bounds are equal, set it at either bound
  end
end

end

