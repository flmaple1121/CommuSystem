# CommuSystem

A MATLAB simulation of a coherent **M-ASK / M-PAM communication link**, developed for a digital communications lab. Explore each transmitter and receiver stage, synchronize frames with a Barker preamble, and compare measured symbol error rates (SER) with the ideal AWGN result.

## Signal chain

```text
Barker-13 preamble + random payload
    -> M-ASK mapping -> upsampling -> RRC pulse shaping
    -> DAC interpolation -> real RF modulation -> AWGN
    -> coherent demodulation -> low-pass filter -> baseband resampling
    -> RRC matched filter -> sampling phase selection -> blind AGC
    -> symbol decisions -> preamble alignment -> payload SER
```

The preamble repeats Barker-13 four times and maps its signs to the outermost constellation levels. The receiver selects a sampling phase by energy, estimates gain from received power, and locates the known preamble by correlation. The carrier frequency and phase are assumed known; carrier recovery and frequency offsets are not modeled.

## Requirements

- MATLAB with **Communications Toolbox** (`rcosdesign`, `awgn`, `finddelay`, constellation and eye diagrams).
- **Signal Processing Toolbox** (`upsample`, `periodogram`, `hann`, and filter-design tools).
- Verified with **MATLAB R2024a**. No Python packages or external downloads are needed.

## Quick start

Clone the repository, open MATLAB, and set its current folder to the repository root:

```sh
git clone https://github.com/flmaple1121/CommuSystem.git
cd CommuSystem
```

Run the complete link with intermediate signal plots, RF power spectral density, an eye diagram, and payload SER:

```matlab
run_lab
```

Run the Monte Carlo sweep and plot simulated versus theoretical SER:

```matlab
run_ser_simulation
```

Both entry points are scripts: edit their configuration sections to change the experiment. Each resets the workspace, closes existing figures, and sets the random seed to 42 for repeatable results.

From a terminal with MATLAB on the path, the equivalent commands are `matlab -batch "run_lab"` and `matlab -batch "run_ser_simulation"`. Use the MATLAB desktop to inspect the figures interactively.

## Default experiments

| Parameter | `run_lab.m` | `run_ser_simulation.m` |
| --- | --- | --- |
| Modulation order | 4-ASK | 16-ASK |
| Symbol rate | 3.84 MBaud | 3.84 MBaud |
| Carrier frequency | 1.92 GHz | 1.92 GHz |
| RF sample rate | 4.5 GHz | 4.5 GHz |
| Baseband samples per symbol | 64 | 64 |
| RRC roll-off / span | 0.22 / 40 symbols | 0.22 / 40 symbols |
| Payload symbols per frame | 2,000 | 2,000 per SNR point |
| Es/N0 | 10 dB | 0 to 35 dB, in 1 dB steps |

Try `M = 2`, `4`, `8`, or `16`. Keep transmitter and receiver parameters consistent. Changing the symbol rate or RF sample rate also requires checking the supplied low-pass filter design. The RF sample rate must cover the carrier and its occupied sidebands.

### Runtime and statistical precision

The RF waveform has approximately **1,172 samples per symbol**. The defaults use short frames to keep memory use practical. Raising `payload_len` to 100,000 creates about 117 million RF samples per frame and can require several gigabytes of memory across processing stages.

With 2,000 payload symbols, one error corresponds to SER = 0.0005. High-SNR points with zero observed errors disappear on a logarithmic plot; they do not establish a zero error probability. Increase the payload length for finer estimates if enough RAM is available. Short frames also make preamble energy and filter transients more significant in the measured-power noise calibration and blind AGC.

### Noise convention and theoretical reference

For the real-valued RF waveform and real noise with two-sided PSD `N0/2`, the sample SNR passed to `awgn` is:

```matlab
sample_snr_db = EsN0_dB - 10*log10(Fs_analog/(2*Rs));
```

The ideal coherent M-PAM / M-ASK symbol error probability is:

```matlab
Ps = 2*(M-1)/M * qfunc(sqrt(6*10.^(EsN0_dB/10)/(M^2-1)));
```

The reference assumes ideal timing and gain. The simulation additionally includes interpolation, finite filters, blind AGC, and preamble synchronization, so finite-frame measurements can differ. Out-of-bounds synchronization in the SER sweep counts as an entire frame of symbol errors.

See MathWorks' [real-signal SNR convention](https://www.mathworks.com/help/comm/ref/comm.awgnchannel-system-object.html) and [PAM error-rate expressions](https://www.mathworks.com/help/comm/ug/analytical-expressions-used-in-berawgn-function-and-bit-error-rate-analysis-app.html). The RF spectrum uses a [windowed periodogram](https://www.mathworks.com/help/signal/ref/periodogram.html), with sample-rate and window-energy normalization, in dB/Hz.

## Repository contents

| File | Purpose |
| --- | --- |
| [`Transmitter.m`](Transmitter.m) | Symbol generation, ASK mapping, rectangular or RRC shaping, DAC interpolation, RF modulation, and transmitter plots |
| [`Receiver.m`](Receiver.m) | Coherent mixing, FIR/IIR low-pass filtering, resampling, RRC matched filtering, blind phase/gain estimation, decisions, and receiver plots |
| [`run_lab.m`](run_lab.m) | Complete 4-ASK lab demonstration, maintained as plain MATLAB source |
| [`run_ser_simulation.m`](run_ser_simulation.m) | 16-ASK SER sweep with Barker preamble alignment and the ideal reference curve |
| [`my_filter_coeffs_FIR.mat`](my_filter_coeffs_FIR.mat) | FIR numerator `Num`, used by both entry points |
| [`my_filter_coeffs.mat`](my_filter_coeffs.mat) | Alternative IIR coefficients `Num` and `Den`; pass both to `Receiver` when experimenting |
| [`lpf.fda`](lpf.fda), [`lpf_butterworth.fda`](lpf_butterworth.fda) | Editable MATLAB Filter Designer sessions |
| [`tests/test_communication_system.m`](tests/test_communication_system.m) | Focused mapping, end-to-end recovery, noise, and PSD regression checks |

The receiver uses an RRC matched filter. Its end-to-end examples therefore use RRC transmit shaping; the transmitter's rectangular mode is intended for separate waveform experiments.

## Validation

Run the focused regression checks from the repository root:

```matlab
results = runtests('tests');
assertSuccess(results);
```

## Local reports

PDF/Word/HTML reports, Live Editor reports (`.mlx`), generated plots, submission ZIPs, and the local duplicate submission folder are excluded by `.gitignore`. The original local files are retained; the Live Editor report is stored in `reports/run_lab.mlx` so it does not shadow the executable `run_lab.m`. A fresh clone contains the code and filter assets required to reproduce the experiments.
