# TRACE — Trilinear Channel Estimation for Tri-Hybrid Architectures

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22923404.svg)](https://doi.org/10.5281/zenodo.22923404)

MATLAB code accompanying:

> A. L. F. de Almeida and M. Liu, *"TRACE: A Trilinear Channel Estimation
> Framework for Tri-Hybrid Architecture,"* IEEE Wireless Communications
> Letters (submitted).

Tri-hybrid beamforming adds an electromagnetic (EM) degree of freedom —
pattern-reconfigurable antennas — on top of RF analog combining and baseband
processing. Cycling the radiation pattern over the pilot dimension turns the
received pilots into a third-order PARAFAC tensor whose three factors are the
EM-domain angular factor, the RF-projected steering factor and the
pilot-weighted fading factor. TRACE fits that tensor by alternating least
squares (ALS) and then recovers the physical channel parameters in closed form.

This package reproduces the three panels of the results figure and is small
enough to read end to end. **It has no dependencies:** base MATLAB, no
toolboxes, no Tensorlab, no Tensor Toolbox. It also runs on GNU Octave.

---

## Quick start

```matlab
setup_trace      % puts core/, baselines/ and reproduce/ on the path
demo_trace       % one realization, full pipeline, ~2 s
```

`demo_trace` prints the per-stage errors for a single channel draw and is the
best entry point for reading the code: every step is annotated with the
equation it implements.

To redraw the paper's figure from the shipped results (seconds, no simulation):

```matlab
make_paper_figure
```

To re-run the simulations from scratch (see *Runtime* below):

```matlab
repro_fig_a_training_dimensions
repro_fig_b_baselines
repro_fig_c_probing_states
make_paper_figure
```

---

## Layout

```
setup_trace.m         path setup and dependency check
demo_trace.m          minimal end-to-end example

core/                 the model and the estimator
  ula_steering.m               a(theta), eq. (2)
  khatrirao.m                  the <> product; read its header on ordering
  trace_unfold.m               Y_(1), Y_(2), Y_(3), eqs. (18)-(20)
  trace_cpdgen.m               Y = [[T_A,T_H,T_S]], eqs. (11),(15)
  trace_design_training.m      F_EMbar and F_RF, eqs. (9),(16)
  trace_generate_clean_data.m  one realization, eqs. (2),(3),(10)-(14)
  trace_add_noise.m            AWGN at a prescribed SNR
  trace_als.m                  ALGORITHM 1, STEPS 1-2, eqs. (17)-(21)
  trace_random_init.m          random start for the multi-start loop
  trace_normalize_factors.m    fixes the CP scaling drift
  trace_prepare_recovery.m     dictionaries for eqs. (23),(24)
  trace_physical_recovery.m    ALGORITHM 1, STEPS 3-7, eqs. (22)-(27)
  trace_channel_nmse.m         the reported metric (Section V)
  trace_channel_energies.m     same, numerator and denominator separately
  trace_eval.m                 all simulation metrics and diagnostics
  trace_best_permutation.m     resolves the CP permutation (diagnostics only)

baselines/            the two comparison methods
  trace_somp_prepare.m         joint EM/angular dictionary (M*G atoms)
  trace_somp_baseline.m        greedy simultaneous OMP
  trace_oracle_baseline.m      oracle-support LS bound

reproduce/            one script per panel, plus the figure renderer
  repro_fig_a_training_dimensions.m
  repro_fig_b_baselines.m
  repro_fig_c_probing_states.m
  make_paper_figure.m

results/              the published data (see results/README.md)
```

The two files to read first are **`core/trace_als.m`** (the tensor fit) and
**`core/trace_physical_recovery.m`** (the parameter recovery). Together they
are Algorithm 1 of the paper.

---

## Notation: paper ↔ code

| Paper | Code | Size | Meaning |
|---|---|---|---|
| `Y` | `data.Y` | N_RF × T × P | pilot tensor, eq. (11) |
| `T_A` | `TA` | N_RF × L | RF-projected steering factor, eq. (12) |
| `T_H` | `TH` | T × L | EM-domain angular factor, eq. (13) |
| `T_S = B` | `TS`, `B_hat` | P × L | pilot-weighted fading factor, eqs. (10),(14) |
| `A(θ)` | `A`, `A_hat` | N_r × L | steering matrix, eq. (2) |
| `H̄_EM` | `HEM`, `HEM_hat` | M × L | one-hot EM signatures, eq. (3) |
| `C` | `alpha`, `alpha_hat` | P × L | stacked path gains, eq. (10) |
| `F̄_EM` | `FEMbar`, `training.FEMbar` | T × M | EM probing matrix, eq. (9) |
| `F_RF` | `FRF`, `training.FRF` | N_RF × N_r | RF combiner |
| `s` | `s` | P × 1 | known pilot sequence |
| `Z_A, Z_H, Z_S` | `ZA, ZH, ZS` | — | ALS regression matrices, eqs. (18)-(20) |
| `⋄` (Khatri–Rao) | `khatrirao` | — | see that file on argument ordering |
| `L = R` | `par.L`, `R` | — | paths = tensor rank (single user) |
| `N_s` | `opts.nstarts` | — | ALS restarts, 5 in the paper |
| `μ(F̄_EM)` | `trainDiag.mu_EM`, `muEM` | — | EM mutual coherence, eq. (16) |

Equation numbers throughout the source refer to the published letter. Every
function header names the equations and Algorithm 1 steps it implements.

---

## What each panel shows

**(a) Training dimensions** — `repro_fig_a_training_dimensions.m`.
Non-compressive setting (N_RF = N_r = 4, T ≥ M = 16, L = 4). Pilot slots beyond
the tensor rank remove the low-SNR threshold effect seen at P = L, and extra EM
probing slots give a uniform gain. At P = L the design rule (29) holds with zero
margin and the estimate still breaks down at low SNR — identifiability is not
conditioning.

**(b) Compressive configuration** — `repro_fig_b_baselines.m`.
N_RF = 3 < N_r = 8 and T = 10 < M = 24, so both reduced domains are genuinely
compressive. The greedy dictionary baseline is more accurate below the crossover
at about 7 dB but saturates at a grid-induced floor; TRACE keeps the SNR slope
of the oracle-support bound and stays within a factor of ~1.4 of it above 10 dB.

**(c) EM probing states** — `repro_fig_c_probing_states.m`.
L = 6, N_RF = 4, P = 8, SNR = 10 dB. Kruskal's condition (29) requires T ≥ 4.
At exactly T = 4 it is satisfied and the estimator still fails: median NMSE 0.82,
38 % of realizations above 1, an order of magnitude off the oracle. The
explanatory variable is the EM coherence of eq. (16) — 0.985 at T = 4 against
0.52 at T = 24. The script controls for codebook luck by nesting: one master
codebook is drawn at T = 24 and smaller designs take its first T rows, with the
physical channel held identical across T within each trial.

All panels plot the **median** over MC = 200 realizations. The median rather
than the mean, because near the low-SNR threshold the distribution is
heavy-tailed and a handful of failed draws (NMSE > 1) would otherwise hide the
behaviour of the typical run. The per-realization matrices are kept in
`results/`, so any other statistic can be recomputed without re-running.

---

## Runtime

On a current laptop, at the published MC = 200:

| Script | Approx. time |
|---|---|
| `demo_trace` | 2 s |
| `repro_fig_a_training_dimensions` | 10–25 min |
| `repro_fig_b_baselines` | 30–60 min |
| `repro_fig_c_probing_states` | 20–40 min |

Set `MC = 20` at the top of any script for a quick look — the curve shapes are
already clear at that point, only the tails move. Panel (b) is the slow one:
its SOMP dictionary has M·G = 24 024 atoms.

---

## Reproducibility, honestly

Each script fixes its seed (`rng(11)`, `rng(22)`, `rng(31)`), so **re-running a
script twice on the same machine gives identical output**.

Re-running will *not* reproduce the shipped `.mat` files bit for bit, for two
reasons, neither of which affects any conclusion:

1. These scripts run ALS only. The original runs also evaluated a Gauss–Newton
   solver (Remark 2 of the paper), which consumed random numbers from the same
   stream. Dropping it shifts the stream, so the random ALS restarts and the
   channel draws from the second realization onward differ.
2. Octave's random number generator differs from MATLAB's for the same seed.

What you get is a statistically equivalent run: at MC = 200 the curves are
visually indistinguishable and every number quoted in Section V reproduces to
within Monte-Carlo error. The exact published data is in `results/`, so
`make_paper_figure` always redraws the figure as it appears in the letter.

---

## Requirements

Base MATLAB R2017b or newer (`vecnorm`). No toolboxes. Three
optional functions are detected at run time and replaced by equivalent
fallbacks when absent — `lsqminnorm` → `pinv`, `matchpairs` → exhaustive
matching, `exportgraphics` → `print`. GNU Octave 6+ works as well.

---

## Citing

If you use this code, please cite the letter (see `CITATION.cff`) and the
archived software:

> de Almeida, A. L. F., & Liu, M. (2026). *TRACE: Trilinear Channel Estimation
> for Tri-Hybrid Architectures (MATLAB)* (v1.0.1) [Computer software]. Zenodo.
> https://doi.org/10.5281/zenodo.22923405

## License

MIT — see `LICENSE`.
