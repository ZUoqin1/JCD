function PlotResults(name, t, trueModes, trueIF, result)
% PLOTRESULTS Display the signal, initialization, and decomposition results.
colors = [214, 39, 40; 23, 86, 179; 0, 158, 115] / 255;
orange = [230, 159, 0] / 255;
K = size(trueModes, 1);
signal = sum(trueModes, 1);
fs = result.settings.fs;
L = result.settings.L;
upperFrequency = 100 * ceil(max(trueIF, [], 'all') / 100);
errors = result.modes - trueModes;
% RE is the error norm divided by the true signal norm, expressed as a percentage.
metrics.IF_RE = 100 * vecnorm(result.IF - trueIF, 2, 2) ./ vecnorm(trueIF, 2, 2);
metrics.Mode_RE = 100 * vecnorm(errors, 2, 2) ./ vecnorm(trueModes, 2, 2);
metrics.Signal_RE = 100 * norm(result.signal - signal) / norm(signal);
fprintf('\n%s\n', name);
disp(table((1:K)', metrics.IF_RE, metrics.Mode_RE, ...
    'VariableNames', {'Mode', 'IF_RE_percent', 'Mode_RE_percent'}));
fprintf('Signal RE: %.6f %%\nDecomposition time: %.3f s\n', metrics.Signal_RE, result.seconds);

folder = fullfile(fileparts(mfilename('fullpath')), 'outputs');
if ~isfolder(folder)
    mkdir(folder);
end
save(fullfile(folder, [name, '.mat']), 't', 'trueModes', 'trueIF', 'result', 'metrics');

f = new_figure([name, ': signal STFT']);
ax = new_axes(f, [.29, .255, .695, .72]);
[frequency, magnitude] = gaussian_stft(signal, fs, L);
valid = frequency >= 0 & frequency <= upperFrequency;
imagesc(ax, [t(1), t(end)], frequency(valid), magnitude(valid, :));
set(ax, 'YDir', 'normal', 'Layer', 'top');
colormap(ax, jet_map());
clim(ax, [0, max(magnitude(valid, :), [], 'all')]);
xlim(ax, [t(1), t(end)]);
ylim(ax, [0, upperFrequency]);
xticks(ax, 0:.2:1);
xlabel(ax, 'Time (s)');
ylabel(ax, 'Frequency (Hz)');
print(f, fullfile(folder, [name, '_STFT.png']), '-dpng', '-r120');

% Use the same first-window spectrum as the IF initializer.
nfft = 65536;
window = 0.5 - 0.5 * cos(2 * pi * (0:L - 1)' / (L - 1));
magnitude = abs(fft(signal(1:L).' .* window, nfft));
magnitude = magnitude(1:nfft / 2 + 1);
magnitude = magnitude / max(magnitude);
frequency = (0:nfft / 2)' * fs / nfft;
f = new_figure([name, ': initial-window spectrum']);
ax = new_axes(f, [.29, .255, .66, .72]);
valid = frequency <= upperFrequency;
plot(ax, frequency(valid), magnitude(valid), 'k', 'LineWidth', 3.9);
for k = 1:K
    % Mark the IF actually used by JCD, in the corresponding mode color.
    peak = result.initialIF(k);
    height = interp1(frequency, magnitude, peak);
    plot(ax, [peak, peak], [0, height], '--', 'Color', colors(k, :), 'LineWidth', 2.75);
    plot(ax, peak, height, 'o', 'Color', 'white', 'MarkerFaceColor', colors(k, :), ...
        'MarkerSize', 15.8, 'LineWidth', 1.75);
end
xlim(ax, [0, upperFrequency]);
ylim(ax, [0, 1.05]);
yticks(ax, 0:.2:1);
ytickformat(ax, '%.1f');
xlabel(ax, 'Frequency (Hz)');
ylabel(ax, 'Magnitude');
print(f, fullfile(folder, [name, '_initialization.png']), '-dpng', '-r120');

f = new_figure([name, ': IF estimates']);
ax = new_axes(f, [.29, .255, .695, .72]);
for k = 1:K
    plot(ax, t, trueIF(k, :), '--', 'Color', colors(k, :), 'LineWidth', 2.8);
    plot(ax, t, result.IF(k, :), 'Color', colors(k, :), 'LineWidth', 3.9);
end
xlim(ax, [t(1), t(end)]);
xticks(ax, 0:.2:1);
xlabel(ax, 'Time (s)');
ylabel(ax, 'Frequency (Hz)');
print(f, fullfile(folder, [name, '_IF.png']), '-dpng', '-r120');

f = new_figure([name, ': mode reconstruction']);
limit = ceil(1.05 * max(abs([trueModes; result.modes; errors]), [], 'all'));
height = .72 / (K + .18 * (K - 1));
for k = 1:K
    bottom = .975 - k * height - (k - 1) * .18 * height;
    ax = new_axes(f, [.29, bottom, .695, height]);
    plot(ax, t, trueModes(k, :), 'k', 'LineWidth', .9);
    plot(ax, t, result.modes(k, :), 'Color', colors(k, :), 'LineWidth', 1.1);
    plot(ax, t, errors(k, :), 'Color', orange, 'LineWidth', 1.4);
    xlim(ax, [t(1), t(end)]);
    ylim(ax, [-limit, limit]);
    xticks(ax, 0:.2:1);
    yticks(ax, [-limit, 0, limit]);
    ylabel(ax, sprintf('M_%d', k));
    if k < K
        ax.XTickLabel = [];
    else
        xlabel(ax, 'Time (s)');
    end
end
print(f, fullfile(folder, [name, '_modes.png']), '-dpng', '-r120');

f = new_figure([name, ': signal reconstruction']);
ax = new_axes(f, [.29, .255, .695, .72]);
plot(ax, t, signal, 'k', 'LineWidth', .9);
plot(ax, t, result.signal, 'Color', [108, 117, 125]/255, 'LineWidth', 1.1);
plot(ax, t, result.signal - signal, 'Color', orange, 'LineWidth', 1.8);
limit = ceil(1.05 * max(abs([signal, result.signal])));
xlim(ax, [t(1), t(end)]);
ylim(ax, [-limit, limit]);
xticks(ax, 0:.2:1);
yticks(ax, [-limit, 0, limit]);
xlabel(ax, 'Time (s)');
ylabel(ax, 'Amplitude');
print(f, fullfile(folder, [name, '_signal.png']), '-dpng', '-r120');
end

function [frequency, magnitude] = gaussian_stft(signal, fs, L)
% Gaussian STFT with one-sample spacing and a 512-point FFT.
analytic = hilbert(signal);
N = numel(signal);
windowLength = 2 * ceil(L / 2);
u = linspace(-1, 1, windowLength);
sigma = .28;
window = (pi * sigma^2)^(-.25) * exp(-u.^2 / (2 * sigma^2));
half = windowLength / 2;
localized = complex(zeros(512, N));
% Each column holds the available windowed samples; unused FFT entries stay zero.
for column = 1:N
    offsets = max(-half, 1 - column):min(half - 1, N - column);
    localized(1:numel(offsets), column) = analytic(column + offsets) .* window(offsets + half + 1);
end
magnitude = abs(fftshift(fft(localized, [], 1), 1));
frequency = (-256:255) * fs / 512;
end

function map = jet_map()
x = linspace(0, 1, 256).';
map = [interp1([0, .35, .66, .89, 1], [0, 0, 1, 1, .5], x), ...
    interp1([0, .125, .375, .64, .91, 1], [0, 0, 1, 1, 0, 0], x), ...
    interp1([0, .11, .34, .65, 1], [.5, 1, 1, 0, 0], x)];
end

function f = new_figure(name)
f = figure('Name', name, 'NumberTitle', 'off', 'Color', 'white', ...
    'Units', 'inches', 'Position', [1, 1, 8.4, 6.8], ...
    'PaperUnits', 'inches', 'PaperPosition', [0, 0, 8.4, 6.8], ...
    'PaperSize', [8.4, 6.8], 'InvertHardcopy', 'off');
end

function ax = new_axes(f, position)
ax = axes(f, 'Units', 'normalized', 'Position', position, ...
    'FontName', 'Arial', 'FontSize', 24.5, 'LineWidth', 2.9, ...
    'Box', 'off', 'TickDir', 'out', 'TickLength', [.014, .014], ...
    'XColor', 'black', 'YColor', 'black');
ax.PositionConstraint = 'innerposition';
hold(ax, 'on');
end
