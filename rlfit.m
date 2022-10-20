%  regParams = RLFIT(r,rl,perc,slLim,glfLim,alfLim,fittingMethod,options)
%
%  DESCRIPTION
%  The function uses constrained non-linear optimisation to minimise the
%  function SUM((RL_PRED - RL).^2), with RL_PRED = SLBES - GLF*log10(R) 
%  - ALF*(R-1). RLFIT will return the three unknowns in RL_PRED, i.e. the best-
%  fit source level (SLBES), geometric loss factor (GLF), and absorption loss 
%  factor (ALF). The three values are returned as fields in output structure
%  REGPARAMS, along with the percentile ranks (PERC), the source levels at those
%  percentile ranks (SLPER), and the standard error of the regression (SLSTD).
% 
%  INPUT VARIABLES
%  - r: vector of ranges [m]
%  - rl: vector of received levels. Undefined sound level metric (e.g. RL can 
%    be in SPLrms, SPL02p, SPLp2p or SEL metrics). Same number of elements as 
%    input range vector R.
%  - perc: vector of percentile ranks over which the regression source levels 
%    in REGPARAMS.SLPER are calculated. Leave empty (PERC = []) for default 
%    values (PERC = [1 5 50 95 99]).
%  - slLim: two-element vector with the bottom and top bounds for the estimated 
%    source level. Leave empty (SLLIM = []) for default values (SLLIM = [0 400])
%  - glfLim: two-element vector with the bottom and top limits for the 
%    estimated geometric (spreading) loss factor. Leave empty (GLFLIM = []) 
%    for default values (GLFLIM = [0 50]).
%  - alfLim: two-element vector with the bottom and top limits for the 
%    estimated absorption loss factor. Leave empty (ALFLIM = []) for default 
%    values (ALFLIM = [0 0.1]).
%  - fittingMethod: character string identifying the optimisation method. The 
%    two methods available use multivariate constrained non-linear approximation to 
%    provide a best-fit solution (in a least squares sense).
%    ¬ 'fminsearchbnd': uses FMINSEARCHBND. No special dependencies. DEFAULT
%    ¬ 'lsqnonlin': uses LSQNONLIN. Requires Optimization Toolbox.
%  - options: structure with information about the optimisation process. 
%    Created with function OPTIMSET(see OPTIMSET for details). FMINSEARCHBND 
%    uses these options: Display, TolX, TolFun, MaxFunEvals, MaxIter, 
%    FunValCheck, PlotFcns, and OutputFcn.
%   
%  OUTPUT VARIABLES
%  - regParams: six-field structure containing the regression parameters.
%    ¬ 'glf': geometric (spreading) loss factor (GLF = 10 for cylindrical 
%       spreading, GLF = 20 for spherical spreading)
%    ¬ 'alf': absorption loss factor
%    ¬ 'slBes': best-fit source level
%    ¬ 'slStd': standard error of regression curve.
%       RL_PRED = SLBES - GLF*log10(R) - ALF*(R-1)
%       SLSTD   = SQRT(SUM((RL - RL_PRED).^2)/LENGTH(RL)))
%    ¬ 'slPer': vector of source levels for the percentile ranks PERC
%    ¬ 'perc': vector of percentile ranges (see input argument PERC).
%
%  FUNCTION CALL
%  regParams = RLFIT(r,rl)
%  regParams = RLFIT(r,rl,perc)
%  regParams = RLFIT(r,rl,perc,slLim)
%  regParams = RLFIT(r,rl,perc,slLim,glfLim)
%  regParams = RLFIT(r,rl,perc,slLim,glfLim,alfLim)
%  regParams = RLFIT(r,rl,perc,slLim,glfLim,alfLim,fittingMethod)
%  regParams = RLFIT(r,rl,perc,slLim,glfLim,alfLim,fittingMethod,options)
%
%  PERC, SLLIM, GLFLIM, ALFLIM, FITTINGMETHOD and OPTIONS can be left empty 
%  ([]). Whenever any of these input arguments are left empty or omitted from 
%  the call, their default values will be used. The list of default values is:
%
%  PERC = [1 5 50 95 99];
%  SLLIM = [0 400]
%  GLFLIM = [0 50]
%  ALFLIM = [0 0.1]
%  FITTINGMETHOD = 'fminsearchbnd'
%  OPTIONS = optimset('Display','off','MaxFunEvals',200,'MaxIter',200,...
%                     'TolFun',1e-4,'TolX',1e-4);
%
%  FUNCTION DEPENDENCIES
%  - PERCENTILE (for FITTINGMETHOD = FMINSEARCHBND)
%  - PERCENTILE, LSQNONLIN (for FITTINGMETHOD = LSQNONLIN)
%
%  TOOLBOX DEPENDENCIES
%  - MATLAB (Core) (for FITTINGMETHOD = FMINSEARCHBND)
%  - MATLAB (Core), Optimization Toolbox (for FITTINGMETHOD = LSQNONLIN)
%
%  EXAMPLE
%  1) Create a received level (RL) curve with +/-5 dB random noise
%  sl = 250;
%  glf = 15;
%  alf = 20e-3;
%  r = 1:1000; % range [m]
%  rln = 10*rand(size(r)) - 5; % received level noise [dB]
%  rl = sl - glf*log10(r) - alf*(r-1) + rln; % received level [dB re 1uPa]
%
%  2) Set the curve fitting parameters
%  slLim = [100 300];
%  glfLim = [5 20];
%  alfLim = [1e-1 1e-4];
%  options = optimset('Display','iter','MaxFunEvals',300,'MaxIter',300,...
%     'TolFun',1e-6,'TolX',1e-6);
%
%  3) Run curve fitting
%  regParams = rlfit(r,rl,[],slLim,glfLim,alfLim,'fminsearchbnd',options);
%  rlPred = regParams.slBes - regParams.glf*log10(r) - regParams.alf*(r-1);
%
%  4) Plot results
%  figure
%  hold on
%  plot(r,rl,'b')
%  plot(r,rlPred,'g','LineWidth',1.5)
%  set(gca,'XScale','log')
%  xlabel('Range [m]')
%  ylabel('RL [dB re 1\muPa]')
%  box on
%  legend('Data','Regression')
%  title(sprintf(['Received Level Curve Fitting \\rm (RL = %0.1f - '...
%      '%0.1flog_{10}(r) - %0.0g(r-1), with RMSE = %0.1f'],regParams.slBes,...
%      regParams.glf,regParams.alf,regParams.slStd))
%
%  See also FMINSEARCHBND, LSQNONLIN, OPTIMSET

%  VERSION 1.0 (01 May 2020)
%  Guillermo Jimenez Arranz
%  email: gjarranz@gmail.com

function regParams = rlfit(r,rl,perc,slLim,glfLim,alfLim,fittingMethod,options)

% Input Arguments
if nargin < 3 || isempty(slLim)
    perc = [1 5 50 95 99];
end
if nargin < 4 || isempty(slLim)
    slLim = [0 400];
end
if nargin < 5 || isempty(glfLim)
    glfLim = [0 50];
end
if nargin < 6 || isempty(alfLim)
    alfLim = [0 0.1];
end
if nargin < 7 || isempty(fittingMethod)
    fittingMethod = 'fminsearchbnd';
end
if nargin < 8 || isempty(options)
    options = optimset('Display','off','MaxFunEvals',200,...
        'MaxIter',200,'TolFun',1e-4,'TolX',1e-4);
end

% Default Values
slLimDefault = [0 400];
glfLimDefault = [5 40];
alfLimDefault = [0 0.05];

% Error Control
if ~isvector(r) || ~isvector(rl) || (length(r) ~= length(rl))
    error('R and RL must be vectors the same size')
end

if ~isempty(slLim)
    if ~isvector(slLim) || ~isnumeric(slLim) || length(slLim) ~= 2 ...
            || any(slLim < 0)
        warning(['SLLIM must be a two-element vector of positive numeric '...
            'values. Default values SLLIM = [%d %d] will be used'],...
            slLimDefault(1),slLimDefault(2))
        slLim = slLimDefault; 
    end
else
    slLim = slLimDefault;
end

if ~isempty(glfLim)
    if ~isvector(glfLim) || ~isnumeric(glfLim) || length(glfLim) ~= 2 || any(glfLim < 0)
        warning(['GLFLIM must be a two-element vector of positive numeric '...
            'values. Default values GLFLIM = [%d %d] will be used'],...
            glfLimDefault(1),glfLimDefault(2))
        glfLim = glfLimDefault; 
    end
else
    glfLim = glfLimDefault;
end

if ~isempty(alfLim)
    if ~isvector(alfLim) || ~isnumeric(alfLim) || length(alfLim) ~= 2 || any(alfLim < 0)
        warning(['ALFLIM must be a two-element vector of positive numeric '...
            'values. Default values ALFLIM = [%d %d] will be used'],...
            alfLimDefault(1),alfLimDefault(2))
        alfLim = alfLimDefault; 
    end
else
    alfLim = alfLimDefault;
end

% Remove NaN, Inf and -Inf values from R and RL
ival = ~isnan(r) & ~isinf(r) & ~isnan(rl) & ~isinf(rl);
r = r(ival);
rl = rl(ival);

% Sort Parameter Limits in Ascending Order
slLim = sort(slLim); % sort source level factor limits in ascending order 
glfLim = sort(glfLim); % sort geometric loss factor limits in ascending order 
alfLim = sort(alfLim); % sort arithmetic loss factor limits in ascending order 
    
% Convert Inputs to Column Vectors
r = r(:); % convert ranges to column vector
rl = rl(:); % received levels to column vector

% Extract Paramter Limits
bMin = alfLim(1);
bMax = alfLim(2);
aMin = glfLim(1);
aMax = glfLim(2);
slMin = slLim(1);
slMax = slLim(2);

switch fittingMethod
    % METHOD 1: Constrained Non-Lin LSQ (Optimisation Toolbox not needed)
    case 'fminsearchbnd' 
        % Curve Fitting
        x0 = [mean(slLim) mean(glfLim) mean(alfLim)]; % starting condition
        lb = [slMin aMin bMin]; % lower bounds
        ub = [slMax aMax bMax]; % upper bounds
        rlfitqe = @(x)sum((rl - (x(1) - x(2)*log10(r) - x(3)*(r-1))).^2);
        regParamsBestFit = fminsearchbnd(rlfitqe,x0,lb,ub,options);
        sl = regParamsBestFit(1);
        glf = regParamsBestFit(2);
        alf = regParamsBestFit(3);

        % Source Level Standard Error
        tlPred = glf.*log10(r) + alf*(r-1); % best-fit regression (TL)
        slPred = rl + tlPred; % source level from all transmission loss points
        slStd = std(slPred); % standard error of best-fit regression 

        % Percentile-Fit Solutions
        nPer = length(perc);
        slPer = nan(size(perc));
        for m = 1:nPer
            slPer(m) = percentile (slPred,perc(m),'near');
        end

        % Output
        regParams.glf = glf;
        regParams.alf = alf;
        regParams.slBes = sl;
        regParams.slStd = slStd;
        regParams.slPer = slPer;
        regParams.perc = perc;
    
    % METHOD 2: Constrained Non-Lin LSQ (Optimisation Toolbox needed)
    case 'lsqnonlin'
        % Curve Fitting
        x0 = [mean(slLim) mean(glfLim) mean(alfLim)]; % starting condition
        lb = [slMin aMin bMin]; % lower bounds
        ub = [slMax aMax bMax]; % upper bounds
        rlfit = @(x) x(1) - x(2)*log10(r) - x(3)*(r-1) - rl;
        x = lsqnonlin(rlfit,x0,lb,ub,options);
        sl = x(1);
        glf = x(2);
        alf = x(3);

        % Source Level Standard Error
        tlPred = glf.*log10(r) + alf*(r-1); % best-fit regression (TL)
        slPred = rl + tlPred; % source level from all transmission loss points
        slStd = std(slPred); % standard error of best-fit regression 
        
        % Find Percentile-Fit Solutions
        nPer = length(perc);
        slPer = nan(size(perc));
        for m = 1:nPer
            slPer(m) = percentile (slPred,perc(m),'near');
        end
        
        % Output
        regParams.glf = glf;
        regParams.alf = alf;
        regParams.slBes = sl;
        regParams.slStd = slStd;
        regParams.slPer = slPer;
        regParams.perc = perc;
end
end

