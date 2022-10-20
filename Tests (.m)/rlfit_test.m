%  1) Create a received level (RL) curve with +/-5 dB random noise
sl = 250;
glf = 15;
alf = 20e-3;
r = 1:1000; % range [m]
rln = 10*rand(size(r)) - 5; % received level noise [dB]
rl = sl - glf*log10(r) - alf*(r-1) + rln; % received level [dB re 1uPa]

%  2) Set the curve fitting parameters
slLim = [100 300];
glfLim = [5 20];
alfLim = [1e-1 1e-4];
options = optimset('Display','iter','MaxFunEvals',300,'MaxIter',300,...
   'TolFun',1e-6,'TolX',1e-6);

%  3) Run curve fitting
regParams = rlfit(r,rl,[],slLim,glfLim,alfLim,'fminsearchbnd',options);
rlPred = regParams.slBes - regParams.glf*log10(r) - regParams.alf*(r-1);

%  4) Plot results
figure
hold on
plot(r,rl,'b')
plot(r,rlPred,'g','LineWidth',1.5)
set(gca,'XScale','log')
xlabel('Range [m]')
ylabel('RL [dB re 1\muPa]')
box on
legend('Data','Regression')
title(sprintf(['Received Level Curve Fitting \\rm (RL = %0.1f - '...
    '%0.1flog_{10}(r) - %0.0g(r-1), with RMSE = %0.1f'],regParams.slBes,...
    regParams.glf,regParams.alf,regParams.slStd))
