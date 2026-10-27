# Jet Chirplet Decomposition

MATLAB demonstrations for the two-mode crossing experiment (Section V-B)
and the three-mode crossing and tangency experiment (Section V-C).

Requires MATLAB and Signal Processing Toolbox.

Set this folder as the MATLAB current folder and run:

```matlab
result1 = Test1;
result2 = Test2;
```

| File | Purpose |
|---|---|
| `Test1.m` | Define the two-mode signal and run JCD |
| `Test2.m` | Define the three-mode signal and run JCD |
| `JCD.m` | Initialize IFs, fit and propagate IF jets, and reconstruct modes |
| `PlotResults.m` | Display results and report relative errors |

Each example displays five figures: the signal STFT, the normalized spectrum
of the initial window, IF estimates, mode waveforms, and signal reconstruction.
Red, blue, and green identify modes in decreasing initial IF order. The
initial-window spectrum is black, with colored markers and dashed lines
indicating the initial IFs. In the IF plot, dashed curves are true IFs and
solid curves are estimates.
In waveform plots, black denotes the true waveform and orange denotes the
reconstruction error; reconstructed modes use their mode colors and the
reconstructed signal is gray.

The `outputs` folder is created when an example runs. It contains the figures
and a MAT file with the signal definitions, numerical results, error metrics,
random seed, and solver settings.

Both examples use a 2 kHz sampling frequency and IF jet order `R = 2`.
`Test1` uses `K = 2`, `L = 101`, and `H = 25`; `Test2` uses `K = 3`,
`L = 175`, and `H = 30`. `L` and `H` are the window length and anchor spacing
in samples. The initial IFs come from the first-window spectrum and the initial
IF derivatives are zero. The initial window has zero proximal weights;
subsequent windows use `[0.1, 0.01, 0]`. Each level allows 24 Gauss-Newton
iterations, with damping `1e-4` and tolerance `1e-10`.

The returned structure contains mode waveforms (`result.modes`), IF estimates
in Hz (`result.IF`), the reconstructed signal (`result.signal`), and fitted
jets at the anchors (`result.jets`). Mode arrays have size `K` by the number
of samples. Each jet contains the IF and its time derivatives in Hz/s^r.
