function result = Test1
% TEST1 Two modes with crossing instantaneous frequencies.
rng(0, 'twister');
fs = 2000;
t = (0:2000) / fs;
amplitude = 1 + 0.5 * cos(2 * pi * t);
% Phase is in cycles, so its time derivative is the IF in Hz.
phase = [0.2 + 532*t - 474*t.^2 + 369*t.^3; ...
         0.8 + 50*t + 525*t.^2 - 300*t.^3];
trueModes = amplitude .* cos(2 * pi * phase);
trueIF = [532 - 948*t + 1107*t.^2; ...
          50 + 1050*t - 900*t.^2];
x = sum(trueModes, 1);

% K: mode count; L: window length; H: anchor spacing; R: IF jet order.
K = 2;
L = 101;
H = 25;
R = 2;
start = tic;
result = JCD(x, fs, K, L, H, R);
result.seconds = toc(start);
result.randomSeed = 0;
PlotResults('Test1', t, trueModes, trueIF, result);
end
