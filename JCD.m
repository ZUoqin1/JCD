function result = JCD(x, fs, K, L, H, R)
% JCD Decompose a signal into K modes using IF jets of order R.
% L is the odd window length; H is the anchor spacing, both in samples.
% Each IF jet contains [f, f', ..., f^(R)], with entry r+1 in Hz/s^r.
% Each row of result.modes and result.IF corresponds to one mode.

validateattributes(x, {'numeric'}, {'real', 'vector', 'finite', 'nonempty'});
validateattributes(fs, {'numeric'}, {'real', 'scalar', 'finite', 'positive'});
validateattributes(L, {'numeric'}, {'scalar', 'integer', 'odd', '>=', 17, '<=', min(numel(x), 65536)});
validateattributes(H, {'numeric'}, {'scalar', 'integer', '>=', 1, '<=', L});
validateattributes(K, {'numeric'}, {'scalar', 'integer', '>=', 1, '<=', L});
validateattributes(R, {'numeric'}, {'scalar', 'integer', '>=', 0, '<=', 3});
x = double(x(:).');
assert(norm(x) > 0, 'JCD:signal', 'The signal must have nonzero energy.');

lambda = [0.1, 0.01, 0];
initialIF = initial_frequencies(x, fs, L, K);
z = hilbert(x);
N = numel(x);
half = (L - 1) / 2;
Th = half / fs;
offset = -half:half;
tau = offset / fs;
window = exp(-0.5 * (offset / (L / 5)).^2).';
% Anchors are window-center sample indices; the last window reaches sample N.
anchors = unique([half + 1:H:N - half, N - half]);
q = R + 1;

% Integrating the IF Taylor series gives phase in cycles; atoms use 2*pi.
basis = zeros(q, L);
IFbasis = zeros(q, L);
for r = 0:R
    basis(r + 1, :) = tau.^(r + 1) / factorial(r + 1);
    IFbasis(r + 1, :) = tau.^r / factorial(r);
end
% Diagonal entries of S_R; dividing jet differences by them normalizes the penalty.
scale = factorial(0:R) ./ Th.^(1:q);
upper = (fs / 2) * factorial(0:R) ./ Th.^(0:R);
lower = -upper;
lower(1) = 0;
upper(1) = fs / 2 - eps(fs / 2) / 2;
jets = zeros(K, q, numel(anchors));
coefficients = complex(zeros(K, numel(anchors)));

% Fit each window, keeping the Taylor prediction fixed during continuation.
for n = 1:numel(anchors)
    y = z(anchors(n) + offset).' .* window;
    if n == 1
        % The initial window has no preceding fitted jet to guide it.
        prediction = [initialIF, zeros(K, R)];
        weights = zeros(size(lambda));
    else
        dt = (anchors(n) - anchors(n - 1)) / fs;
        % Taylor transport advances each derivative to the new window center.
        prediction = jets(:, :, n - 1);
        for r = 0:R
            for s = r + 1:R
                prediction(:, r + 1) = prediction(:, r + 1) + ...
                    jets(:, s + 1, n - 1) * dt^(s - r) / factorial(s - r);
            end
        end
        weights = lambda;
    end
    current = min(max(prediction, lower), upper);
    % Each later level starts from the preceding level's fit in this same window.
    for weight = weights
        [current, c] = fit_window(y, window, basis, current, prediction, ...
            scale, Th, lower, upper, weight);
    end
    jets(:, :, n) = current;
    coefficients(:, n) = c;
end

% Overlap-add uses squared window weights for both waveforms and IFs.
modes = zeros(K, N);
IF = zeros(K, N);
overlap = zeros(1, N);
weight = window.'.^2;
for n = 1:numel(anchors)
    index = anchors(n) + offset;
    localModes = coefficients(:, n) .* exp(2i * pi * (jets(:, :, n) * basis));
    modes(:, index) = modes(:, index) + real(localModes) .* weight;
    IF(:, index) = IF(:, index) + (jets(:, :, n) * IFbasis) .* weight;
    overlap(index) = overlap(index) + weight;
end
modes = modes ./ overlap;
IF = IF ./ overlap;

% Display modes in decreasing initial IF order; their indices stay fixed.
order = K:-1:1;
result.modes = modes(order, :);
result.IF = IF(order, :);
result.signal = sum(modes, 1);
result.jets = jets(order, :, :);
result.coefficients = coefficients(order, :);
result.anchors = anchors;
result.initialIF = initialIF(order);
result.settings = struct('fs', fs, 'K', K, 'L', L, 'H', H, 'R', R, ...
    'lambda', lambda, 'initialWindowLambda', [0, 0, 0], ...
    'damping', 1e-4, 'tolerance', 1e-10, 'maxIterations', 24);
end

function initialIF = initial_frequencies(x, fs, L, K)
% Select the K largest eligible spectral peaks in the first window.
nfft = 65536;
% Hann weighting and zero padding give the spectrum used for initialization.
window = 0.5 - 0.5 * cos(2 * pi * (0:L - 1)' / (L - 1));
magnitude = abs(fft(x(1:L).' .* window, nfft));
magnitude = magnitude(1:nfft / 2 + 1);
frequency = (0:nfft / 2)' * fs / nfft;
% Peak spacing is about 2*fs/L Hz; prominence is relative to the largest magnitude.
[~, peaks] = findpeaks(magnitude, ...
    'MinPeakDistance', max(1, floor(2 * nfft / L)), ...
    'MinPeakProminence', 0.08 * max(magnitude));
peaks = peaks(frequency(peaks) > 20 & frequency(peaks) < fs / 2 - 20);
assert(numel(peaks) >= K, 'JCD:initialization', 'The first window must contain at least K eligible peaks.');
[~, order] = sort(magnitude(peaks), 'descend');
initialIF = sort(frequency(peaks(order(1:K))));
end

function [jets, c] = fit_window(y, window, basis, jets, prediction, S, Th, lo, hi, lambda)
% Variable projection eliminates complex coefficients at every trial jet.
[K, q] = size(jets);
P = K * q;
lower = repmat(lo(:), K, 1);
upper = repmat(hi(:), K, 1);
% Pack all coordinates of mode 1, then mode 2, and so on.
theta = min(max(reshape(jets.', [], 1), lower + eps), upper - eps);
reference = reshape(prediction.', [], 1);
scale = repmat(S(:), K, 1);
% The penalty divides deviations by S_R; numerical steps use Th*S_R.
numericalScale = Th * scale;
penaltyScale = sqrt(lambda) ./ scale;
ynorm = max(norm(y), eps);
eta = 1e-4;
tolerance = 1e-10;
[residual, c, A, e, gram] = evaluate(theta);
objective = 0.5 * (residual.' * residual);

for iteration = 1:24
    % Exact Jacobian of the reduced waveform residual, followed by the penalty.
    Ai = (gram \ eye(K)) * A';
    Jc = complex(zeros(numel(y), P));
    for k = 1:K
        for r = 1:q
            % Coordinate r represents IF derivative order r-1.
            derivative = (2i * pi * basis(r, :).') .* A(:, k);
            v = derivative * c(k);
            % Remove changes absorbed by coefficient fitting, then include
            % the residual-dependent term in the variable-projection Jacobian.
            projected = v - A * (Ai * v);
            adjoint = zeros(K, 1);
            adjoint(k) = derivative' * e;
            Jc(:, (k - 1) * q + r) = -(projected + Ai' * adjoint) / ynorm;
        end
    end
    % Appending the penalty Jacobian adds lambda*W to the normal-equation matrix.
    J = [real(Jc); imag(Jc); diag(penaltyScale)];
    J = J .* numericalScale.';
    u = -((J.' * J + eta * eye(P)) \ (J.' * residual));
    if any(~isfinite(u)) || norm(u) <= tolerance
        break
    end
    % Convert the normalized step u back to the original jet coordinates.
    increment = numericalScale .* u;
    accepted = false;
    % Refit coefficients at each trial; accept the first sufficient objective decrease.
    for gamma = [1, 1/2, 1/4, 1/8, 1/16]
        trial = min(max(theta + gamma * increment, lower), upper);
        [trialResidual, trialC, trialA, trialE, trialGram] = evaluate(trial);
        trialObjective = 0.5 * (trialResidual.' * trialResidual);
        if isfinite(trialObjective) && objective - trialObjective >= tolerance
            theta = trial;
            residual = trialResidual;
            objective = trialObjective;
            c = trialC;
            A = trialA;
            e = trialE;
            gram = trialGram;
            accepted = true;
            break
        end
    end
    if ~accepted
        break
    end
end
jets = reshape(theta, q, K).';

    function [residual, c, A, e, gram] = evaluate(parameters)
        localJets = reshape(parameters, q, K).';
        % .': transpose without conjugation; A' below is the conjugate transpose.
        A = exp(2i * pi * (localJets * basis)).' .* window;
        gram = A' * A;
        % Complex coefficients represent amplitude and phase offset at the anchor.
        c = gram \ (A' * y);
        e = y - A * c;
        % Half the squared residual norm equals the normalized waveform-fitting objective
        % plus the proximal penalty; the final level has lambda = 0.
        residual = [real(e / ynorm); imag(e / ynorm); ...
            penaltyScale .* (parameters - reference)];
    end
end
