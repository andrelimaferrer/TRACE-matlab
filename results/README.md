# Published results

The data behind the three panels of the results figure. These files are what
`make_paper_figure` reads, so the figure can be redrawn without re-running any
simulation.

Each file keeps the **per-realization** NMSE matrices (rows = operating point,
columns = the MC = 200 Monte-Carlo draws), not just the averaged curves. That
is deliberate: the paper reports medians and NMSE > 1 failure rates, and
neither can be recovered from a mean curve. Every file also carries a
`solver` and a `provenance` string.

| File | Panel | Key variables |
|---|---|---|
| `figure_a_data.mat` | (a) training dimensions | `SNRdB`, `Cfg`, `labels`, `all_nmse{c}` (8 × 200 each), `curves` |
| `figure_b_data.mat` | (b) compressive configuration | `SNRdB`, `als`, `somp`, `oracle` (10 × 200 each), `med_*`, `par` |
| `figure_c_data.mat` | (c) EM probing states | `Tvec`, `muEM`, `ratALS`/`ratSOMP`/`ratOracle` (6 × 200), `numALS`/`denTrue`, `med*`, `Tmin` |

The original run files also stored the full MATLAB workspace (one was 21 MB);
these are stripped to the arrays that matter, about 60 KB each.

## Note on panel (a)

`figure_a_data.mat` is the published panel-(a) data, produced with the
Gauss–Newton solver at N_s = 2 restarts (`solver` field records this). The
letter's Remark 2 states the two solvers were indistinguishable at these
operating points, and `repro_fig_a_training_dimensions.m` runs the ALS path at
N_s = 5 like the other two panels. Re-run that script to regenerate this file
with the ALS solver; the numbers move by less than the Monte-Carlo error.

## Recomputing the quoted numbers

```matlab
F = load('results/figure_c_data.mat');
t = @(v) find(F.Tvec==v,1);
median(F.ratALS(t(4),:))              % 0.82
100*mean(F.ratALS(t(4),:) > 1)        % 38 %
F.muEM(t(4)) / F.muEM(t(24))          % 0.985 / 0.52
```
