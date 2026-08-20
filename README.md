# DO-2025: Simulation Code for "From Phase Ib/II to Seamless Phase II/III: A Simulation-Based Discussion of Design Strategies for Dose Optimization and Confirmatory Phase"

This repository contains the R simulation code used to generate every figure and table in the manuscript and Supplementary Materials.

## Repository structure

| Folder | Strategy |
|---|---|
| `phase 1b2/` | Strategy 1 — DO in Phase Ib/II design |
| `phase 23 op/` | Strategy 2 — DO in operational seamless Phase II/III design |
| `phase 23 op calibration/` | Calibration scripts for the operational seamless design's decision thresholds |
| `phase 23 inf/` | Strategy 3 — DO in inferential seamless Phase II/III design |

Each folder contains the scenario configuration, simulation driver, and post-processing scripts for that strategy.

## Reproducing a scenario

Each scenario in Table 1 (dose-response) and the OS0/OS1 survival scenarios is set by editing the following parameters at the top of the simulation driver script in each folder:

```r
# Dose-response (efficacy/toxicity) parameters
p_ORR_h = <high-dose ORR>
p_TOX_h = <high-dose TOX>
p_ORR_l = <low-dose ORR>
p_TOX_l = <low-dose TOX>
p_ORR_c = <control ORR>
p_TOX_c = <control TOX>

# Survival parameters
lambda_c = 0.06          # control hazard rate (median OS = 11.6 months)
HR = c(HR_low, HR_high)  # hazard ratios for low- and high-dose arms vs. control
```

### Table 1 scenarios (dose-response)

| Scenario | Control ORR | Control TOX | LD ORR | LD TOX | HD ORR | HD TOX | Description |
|---|---|---|---|---|---|---|---|
| S01 | 0.10 | 0.10 | 0.10 | 0.10 | 0.10 | 0.10 | LD and HD equally ineffective but safe |
| S02 | 0.10 | 0.10 | 0.10 | 0.30 | 0.10 | 0.30 | LD and HD equally ineffective and unsafe |
| S1  | 0.10 | 0.10 | 0.25 | 0.10 | 0.30 | 0.10 | HD more effective; both safe |
| S2  | 0.10 | 0.10 | 0.25 | 0.05 | 0.30 | 0.10 | HD more effective; LD safer |
| S3  | 0.10 | 0.10 | 0.25 | 0.10 | 0.25 | 0.10 | LD and HD equally effective and safe |
| S4  | 0.10 | 0.10 | 0.25 | 0.05 | 0.25 | 0.10 | LD and HD equally effective; LD safer |

Map each row to the script parameters as: `p_ORR_c` = Control ORR, `p_TOX_c` = Control TOX, `p_ORR_l` = LD ORR, `p_TOX_l` = LD TOX, `p_ORR_h` = HD ORR, `p_TOX_h` = HD TOX.

S01–S02 are used for Type I error evaluation (paired with OS0); S1–S4 are used for power evaluation (paired with OS1).

### Survival scenarios

| Scenario | Control median OS | HR (low dose) | HR (high dose) | Description |
|---|---|---|---|---|
| OS0 | 11.6 months (`lambda_c = 0.06`) | 1.0 | 1.0 | No survival benefit for either dose |
| OS1 | 11.6 months (`lambda_c = 0.06`) | 0.70 | 0.70 | Both doses improve median OS to 16.6 months |

Additional survival scenarios used in Appendix D are:

| Scenario | Control median OS | HR (low dose) | HR (high dose) | Description |
|---|---|---|---|---|
| OS2 | 11.6 months (`lambda_c = 0.06`) | 0.80 | 0.60 | High dose has more survival benefit than low dose |
| OS3 | 11.6 months (`lambda_c = 0.06`) | 1.0 | 0.70 | Only high dose has survival benefit |

