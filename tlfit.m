%  regParams = TLFIT(r,tl,perc,glfLim,alfLim,fittingMethod,options)
%
%  DESCRIPTION
%  The function uses constrained non-linear optimisation to minimise the
%  function SUM((TL_PRED - TL).^2), with TL_PRED = GLF*log10(R) + ALF*(R-1). 
%  TLFIT will return the two unknowns in TL_PRED, i.e. the geometric loss 
%  factor (GLF) and absorption loss factor (ALF). The two values are returned 
%  as fields in output structure REGPARAMS, along with the standard error of 
%  the regression (TLSTD).
% 
%  INPUT VARIABLES
%  - r: vector of ranges [m]
%  - tl: vector of transmission losses [dB]. Same number of elements as input 
%    range vector R.
%  - glfLim: two-element vector with the bottom and top limits for the 
%    estimated geometric (spreading) loss factor. Leave empty (GLFLIM = []) 
%    for default values (GLFLIM = [0 50]).
%  - alfLim: two-element vector with the bottom and top limits for the 
%    estimated absorption loss factor. Leave empty (ALFLIM = []) for default 
%    values (ALFLIM = [0 0.1]).
%  - fittingMethod: character string identifying the optimisation method. The 
%    two methods available use multivariate constrained non-linear approximation
%    to provide a best-fit solution (in a least squares sense).
%    ¬ 'fminsearchbnd': uses FMINSEARCHBND. No special dependencies. DEFAULT
%    ¬ 'lsqnonlin': uses LSQNONLIN. Requires Optimization Toolbox.
%  - options: structure with information about the optimisation process. 
%    Created with function OPTIMSET(see OPTIMSET for details). FMINSEARCHBND 
%    uses these options: Display, TolX, TolFun, MaxFunEvals, MaxIter, 
%    FunValCheck, PlotFcns, and OutputFcn.
%   
%  OUTPUT VARIABLES
%  - regParams: three-field structure containing the regression parameters.
%    ¬ 'glf': geometric (spreading) loss factor (GLF = 10 for cylindrical 
%       spreading, GLF = 20 for spherical spreading)
%    ¬ 'alf': absorption loss factor
%    ¬ 'tlStd': standard error of regression curve.
%       TL_PRED = GLF*log10(R) + ALF*(R-1)
%       SLSTD   = SQRT(SUM((TL - TL_PRED).^2)/LENGTH(TL)))
%
%  FUNCTION CALL
%  regParams = TLFIT(r,tl,glfLim,alfLim,fittingMethod,options)
%  regParams = TLFIT(r,tl)
%  regParams = TLFIT(r,tl,glfLim)
%  regParams = TLFIT(r,tl,glfLim,alfLim)
%  regParams = TLFIT(r,tl,glfLim,alfLim,fittingMethod)
%  regParams = TLFIT(r,tl,glfLim,alfLim,fittingMethod,options)
%
%  GLFLIM, ALFLIM, FITTINGMETHOD and OPTIONS can be left empty ([]). Whenever 
%  any of these input arguments are left empty or omitted from the call, their 
%  default values will be used. The list of default values is:
%
%  GLFLIM = [0 50]
%  ALFLIM = [0 0.1]
%  FITTINGMETHOD = 'fminsearchbnd'
%  OPTIONS = optimset('Display','off','MaxFunEvals',200,'MaxIter',200,...
%                     'TolFun',1e-4,'TolX',1e-4);
%
%  FUNCTION DEPENDENCIES
%  - None (for FITTINGMETHOD = FMINSEARCHBND)
%  - LSQNONLIN (for FITTINGMETHOD = LSQNONLIN)
%
%  TOOLBOX DEPENDENCIES
%  - MATLAB (Core) (for FITTINGMETHOD = FMINSEARCHBND)
%  - MATLAB (Core), Optimization Toolbox (for FITTINGMETHOD = LSQNONLIN)
%
%  EXAMPLE
%  1) Create a transmission loss (TL) curve with +/-5 dB random noise
%  glf = 15;
%  alf = 20e-3;
%  r_dec = [1 25:25:1000]; % range [m]
%  tln = 5*rand(size(r_dec)) - 2.5; % transmission loss (noise) [dB]
%  tls = glf*log10(r_dec) + alf*(r_dec-1); % transmission loss (clean) [dB]
%  tl_dec = tls + tln; % received level [dB re 1uPa]
%  r = 1:1000;
%  tl = interp1(r_dec,tl_dec,r,'pchip');
%
%  2) Set the curve fitting parameters
%  glfLim = [5 20];
%  alfLim = [1e-1 1e-4];
%  options = optimset('Display','iter','MaxFunEvals',300,'MaxIter',300,...
%     'TolFun',1e-6,'TolX',1e-6);
%
%  3) Run curve fitting
%  regParams = tlfit(r,tl,glfLim,alfLim,'fminsearchbnd',options);
%  tlPred = regParams.glf*log10(r) + regParams.alf*(r-1);
%
%  4) Plot results
%  figure
%  hold on
%  plot(r,tl,'b')
%  plot(r,tlPred,'g','LineWidth',1.5)
%  set(gca,'XScale','log')
%  xlabel('Range [m]')
%  ylabel('TL [dB]')
%  box on
%  legend('Data','Regression')
%  title(sprintf(['Transmission Loss Curve Fitting \\rm (TL = '...
%      '%0.1flog_{10}(r) - %0.0g(r-1), with RMSE = %0.1f'],...
%      regParams.glf,regParams.alf,regParams.tlStd))
%
%  See also FMINSEARCHBND, LSQNONLIN, OPTIMSET

%  VERSION 1.0 (01 May 2020)
%  Guillermo Jimenez Arranz
%  email: gjarranz@gmail.com

function regParams = tlfit(r,tl,glfLim,alfLim,fittingMethod,options)

% Input Arguments
if nargin < 3 || isempty(glfLim)
    glfLim = [0 50];
end
if nargin < 4 || isempty(alfLim)
    alfLim = [0 0.1];
end
if nargin < 5 || isempty(fittingMethod)
    fittingMethod = 'fminsearchbnd';
end
if nargin < 6 || isempty(options)
    options = optimset('Display','off','MaxFunEvals',200,...
        'MaxIter',200,'TolFun',1e-4,'TolX',1e-4);
end

% Default Values
glfLimDefault = [5 40];
alfLimDefault = [0 0.05];

% Error Control
if ~isvector(r) || ~isvector(tl) || (length(r) ~= length(tl))
    error('R and RL must be vectors the same size')
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

% Remove NaN, Inf and -Inf values from R and TL
ival = ~isnan(r) & ~isinf(r) & ~isnan(tl) & ~isinf(tl);
r = r(ival);
tl = tl(ival);

% Sort Parameter Limits in Ascending Order
glfLim = sort(glfLim); % sort geometric loss factor limits in ascending order 
alfLim = sort(alfLim); % sort arithmetic loss factor limits in ascending order 
    
% Convert Inputs to Column Vectors
r = r(:); % convert ranges to column vector
tl = tl(:); % transmission loss to column vector

% Extract Paramter Limits
bMin = alfLim(1);
bMax = alfLim(2);
aMin = glfLim(1);
aMax = glfLim(2);

switch fittingMethod
    % METHOD 1: Constrained Non-Lin LSQ (Optimisation Toolbox not needed)
    case 'fminsearchbnd' 
        % Curve Fitting
        x0 = [mean(glfLim) mean(alfLim)]; % starting condition
        lb = [aMin bMin]; % lower bounds
        ub = [aMax bMax]; % upper bounds
        tlfitqe = @(x)sum((tl - (x(1)*log10(r) + x(2)*(r-1))).^2);
        regParamsBestFit = fminsearchbnd(tlfitqe,x0,lb,ub,options);
        glf = regParamsBestFit(1);
        alf = regParamsBestFit(2);

        % Source Level Standard Error
        tlPred = glf.*log10(r) + alf*(r-1); % best-fit regression (TL)
        tlStd = sqrt(sum((tl - tlPred).^2)/length(tl)); % std error of regression 

        % Output
        regParams.glf = glf;
        regParams.alf = alf;
        regParams.tlStd = tlStd;
    
    % METHOD 2: Constrained Non-Lin LSQ (Optimisation Toolbox needed)
    case 'lsqnonlin'
        % Curve Fitting
        x0 = [mean(glfLim) mean(alfLim)]; % starting condition
        lb = [aMin bMin]; % lower bounds
        ub = [aMax bMax]; % upper bounds
        tlfit = @(x) x(1)*log10(r) - x(2)*(r-1) - tl;
        x = lsqnonlin(tlfit,x0,lb,ub,options);
        glf = x(1);
        alf = x(2);

        % Source Level Standard Error
        tlPred = glf.*log10(r) + alf*(r-1); % best-fit regression (TL)
        tlStd = sqrt(sum((tl - tlPred).^2)/length(tl)); % std error of regression 
        
        % Output
        regParams.glf = glf;
        regParams.alf = alf;
        regParams.tlStd = tlStd;
end
end

