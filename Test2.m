function result = Test2
% TEST2 Three modes with crossing and tangent instantaneous frequencies.
rng(0, 'twister');
fs = 2000;
t = (0:2000) / fs;
% Rows correspond to the decreasing, sinusoidal, and constant IFs, respectively.
amplitude = [0.5 + 0.2*cos(2*pi*t + pi/4); ...
             1 + 0.10*cos(2*pi*t + 3*pi/4); ...
             1 + 0.15*cos(2*pi*t + 7*pi/4)];
% Phase is in cycles, so its time derivative is the IF in Hz.
phase = [(400 + 80*cos(3*pi/20))*t + ...
             200*pi*sin(3*pi/20)*(43*t/50 - t.^2) + 5/16; ...
         400*t + (16/pi)*sin(5*pi*t); ...
         320*t + 1/4];
trueModes = amplitude .* cos(2 * pi * phase);
trueIF = [400 + 80*cos(3*pi/20) + 200*pi*sin(3*pi/20)*(43/50 - 2*t); ...
          400 + 80*cos(5*pi*t); ...
          320*ones(size(t))];
x = sum(trueModes, 1);

% K: mode count; L: window length; H: anchor spacing; R: IF jet order.
K = 3;
L = 175;
H = 30;
R = 2;
start = tic;
result = JCD(x, fs, K, L, H, R);
result.seconds = toc(start);
result.randomSeed = 0;
PlotResults('Test2', t, trueModes, trueIF, result);
end
