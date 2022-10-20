% Create a transmission loss (TL) curve with +/-5 dB random noise
glf = 15;
alf = 20e-3;
r_dec = [1 25:25:1000]; % range [m]
tln = 5*rand(size(r_dec)) - 2.5; % transmission loss (noise) [dB]
tls = glf*log10(r_dec) + alf*(r_dec-1); % transmission loss (clean) [dB]
tl_dec = tls + tln; % received level [dB re 1uPa]
r = 1:1000;
tl = interp1(r_dec,tl_dec,r,'pchip');

% Set the curve fitting parameters
glfLim = [5 20];
alfLim = [1e-1 1e-4];
options = optimset('Display','iter','MaxFunEvals',300,'MaxIter',300,...
   'TolFun',1e-6,'TolX',1e-6);

% Run curve fitting
regParams = tlfit(r,tl,glfLim,alfLim,'fminsearchbnd',options);
tlPred = regParams.glf*log10(r) + regParams.alf*(r-1);

% Plot results
figure
hold on
plot(r,tl,'b')
plot(r,tlPred,'g','LineWidth',1.5)
set(gca,'XScale','log')
xlabel('Range [m]')
ylabel('TL [dB]')
box on
legend('Data','Regression')
title(sprintf(['Transmission Loss Curve Fitting \\rm (TL = '...
    '%0.1flog_{10}(r) - %0.0g(r-1), with RMSE = %0.1f'],...
    regParams.glf,regParams.alf,regParams.tlStd))