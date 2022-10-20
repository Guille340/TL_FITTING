
function sl = percfit(rtl,tl,rrl,rl,perc)

% DESCRIPTION: Best Fit of Simulated RL (i.e. SL - TL) to Measured RL

tli = interp1(rtl,tl,rrl);
sl = rl + tli;
M = length(perc);
slp = zeros(1,M);
for m = 1:M
    slp(m) = percentile (sl,perc(m),'near');
end
sl = slp;




