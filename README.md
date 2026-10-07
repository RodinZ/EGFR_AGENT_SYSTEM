# High-Affinity, pH-Selective EGFR Binder Design: Technical Route and Screening Report

---

## 1. Objectives and Design Principles

The goal is to design high-affinity binders against EGFR, ranked by the following priorities:

| Priority | Objective | Criterion |
|---|---|---|
| **P1** | **pH selectivity** | binds human EGFR at pH 6.5, no detectable binding at pH 7.4 |
| **P2** | **mouse cross-reactivity** | the same sequence also binds mouse EGFR, enabling later mouse-model validation |
| **P3** | **human EGFR affinity** | the stronger the better, ranked below the two above |

---

## 2. Technical Route

![Workflow for de novo EGFR binder design](fig_agent_workflow_en.svg)

*Figure 1　Two-stage workflow. Stage 1 derives hotspot constraints from the target structure and generates de novo backbones along two parallel routes; Stage 2 applies multi-model structural filtering and three-species KD scoring, then pH post-design, emitting two candidate pools.*

The route is a closed loop of epitope analysis → guidance set and SFT → generation → structural filtering → KD scoring → pH post-design → panel selection. Human/mouse conserved epitopes serve as hotspot constraints, and a guidance set built from experimental measurements plus the agent's own scored candidates supervises a multi-objective fine-tune of the generator; RFdiffusion3 and BindCraft run in parallel to produce ~90 aa mini-binders; candidates pass multi-model structural filtering and ESM2-3B three-species scoring before being tiered into the pool; His centres are then placed at the interface in pH post-design, and every new sequence is re-scored for structure and cross-species probability. The loop converges on candidates combining pH selectivity, mouse cross-reactivity and high affinity.

Execution units, referenced by code below:

| Code | Stage | Section |
|---|---|---|
| A1 | epitope analysis | §3 |
| A5 | guidance-set scoring | §4.1 |
| A2′ | multi-objective end-to-end SFT | §4.2 |
| A2 | dual-route generation | §5.1 |
| A3 | candidate-pool maintenance | §5.3 |
| A4 | multi-model structural filtering | §6.1 |
| A5 | three-species KD scoring | §6.2 |
| A6 | pH post-design | §6.3 |
| A7 | panel selection and submission | §7.2 / §8 |

Two workflow-level constraints: the generating units (A2, A6) are separate from the rejecting units (A4, A5), and A4's thresholds are neither visible nor adjustable to A2; low-confidence candidates are demoted rather than deleted, preserving traceable re-entry.

<details>
<summary>Text version of the flow diagram</summary>

```
  EGFR target data (human / mouse / cyno)
        ↓
  A1 epitope analysis: conserved + divergent epitopes → hotspot constraints
        ↓
  A5  guidance set: experimental data + agent-generated candidates
                    → three-species scoring → weighted merged set
        ↓
  A2′ multi-objective SFT: objective-conditioned fine-tune (P1 > P2 > P3)
                    → novelty filter against the guidance set
        ↓
  A2 generation ──┬── Route A  RFdiffusion3 → SolubleMPNN → 9 rounds partial diffusion
                  └── Route B  BindCraft / FreeBindCraft → SolubleMPNN → refold
        ↓
  A4 multi-model structural filtering  i_pTM / pTM / i_pAE / Unbound_pLDDT
        ↓
  A5 KD predictor  ESM2-3B + three-species heads
        ↓   ⟲ A6 output re-scored through A4 / A5
  A6 pH post-design: interface His centres + pH 6.5 composite ranking
        ↓
  A7 panel selection: global_top + diverse_panel → submission CSV
```

</details>

---

## 3. Design Basis: EGFR Epitope Analysis

Contact epitopes between the cetuximab-like antibody and human EGFR were defined from the 1YY9 structure (EGFR chain A; Fab light / heavy chains C/D; 4.5 Å atomic contact cutoff), then mapped onto mouse EGFR to separate conserved from divergent epitopes:

| Epitope type | Positions | Design implication |
|---|---|---|
| conserved (shared human / mouse) | 373, 406, 408, 432, 433, 435, 436, 439, 441, 462, 464, 465, 489, 490, 493 | fixed as anchors, preserving the original binding mode |
| divergent (human / mouse differ) | 377 R/K, 442 S/G, 467 K/R, 491 I/M, 492 S/N, 495 G/A, 497 N/K | side chains take small or polar residues compatible with both species |

Sources: PDB 1YY9 + UniProt human EGFR P00533 / mouse EGFR Q01279 | Criterion: 4.5 Å atomic contact cutoff

**Numbering**: positions above are UniProt P00533 precursor numbering, which includes the 24-residue signal peptide. PDB 1YY9 chain A uses mature numbering, so subtract 24 before indexing into the structure — conserved 373/406/408 map to 349/382/384, divergent 377/467/492 map to 353/443/468. Hotspot indices handed to a generator must be in the numbering of the structure file it reads; a 24-residue offset silently moves the hotspot to an unrelated surface patch.

Hotspots are drawn mainly from conserved epitopes, so that P2 (cross-species binding) is constrained structurally at generation time rather than filtered for afterwards.

---

## 4. Guidance Set and Multi-Objective SFT

Generation is conditioned on a guidance set built from two sources and consumed by an end-to-end fine-tune, rather than on hand-written heuristics alone.

### 4.1 Guidance Set Construction

| Source | Content | Label |
|---|---|---|
| Experimental | ProteinBase EGFR 833 + Claude 90 (human), mouse 90, cyno 90 | measured binder / non-binder, KD where available |
| Agent-generated | the 492 de novo candidates of §5.2 | structural metrics plus ESM2-3B three-species scores, used as pseudo-labels |

Both sources are re-scored with the ESM2-3B three-species heads and merged into one table carrying sequence, source, label or pseudo-label, per-species probability, structural metrics where present, and a confidence weight.

The experimental set alone is small — 923 human records and under 100 each for mouse and cyno — and skewed towards antibody-like binders. The agent-generated set covers the de novo sequence distribution the generator actually samples from, which the experimental set does not. Experimental records therefore carry a higher confidence weight than pseudo-labels, so the latter broaden coverage without dominating the objective.

Outputs: `gene_data/guidance_set.csv` (merged, scored, weighted) and `gene_data/guidance_priors.json` (position-specific residue bias, length and net-charge windows, liability pre-filter thresholds).

### 4.2 Multi-Objective End-to-End SFT

The guidance set supervises an end-to-end fine-tune of the generator, so the three objectives steer sampling directly instead of being applied only as a post-hoc filter.

- **Conditioning**: objective embeddings for pH selectivity, cross-species binding and affinity, so a design can be sampled under a chosen objective mix
- **Loss**: sequence likelihood plus three weighted objective terms, ordered P1 > P2 > P3 to match §1
- **Sample weighting**: experimental records above agent pseudo-labels, as in §4.1
- **Output**: a fine-tuned generator used by A2 alongside unconditioned sampling, with both kept in the pool so the fine-tune cannot silently narrow the search

**Novelty filter**: every sampled design is checked for sequence identity against the full guidance set, and anything above the identity ceiling is discarded. Training on public measurements is the same practice that produced the KD predictor, but a fine-tuned model can memorise and re-emit training sequences; the filter is what keeps that from turning into an existing binder reappearing as a design.

**Compliance boundary**: only distribution-level statistics, weights and model parameters are passed downstream. No sequence, fragment or motif from an existing binder may seed generation, and no sequence in the guidance set enters the candidate pool or the submission.

---

## 5. Candidate Generation

### 5.1 De Novo Mini-Binder Generation (Dual Route)

Roughly 90 aa mini-binders are generated de novo against EGFR along two complementary routes, run in parallel without filtering each other:

| Route | Pipeline | Current output |
|---|---|---|
| **Route A: RFdiffusion3** | RFdiffusion3 backbone → SolubleMPNN → 9 rounds partial diffusion → multi-model filtering | quick screen, 32 designs (no structural metrics yet) |
| **Route B: BindCraft** | FreeBindCraft / BindCraft → SolubleMPNN or co-design → refold / structural filtering / ranked | 87 structure-pass, 32 rank candidates |

Inverse folding uses SolubleMPNN throughout.

### 5.2 Candidate Summary

Candidates are collected in `gene_data/biogen.csv`: **492 records (373 unique sequences)**, tiered by how directly they can be acted on:

| Tier | Rank range | Count |
|---|---|---|
| `structure_pass_best` | 1–87 | 87 |
| `bindcraft_ranked` | 88–119 | 32 |
| `rf_sequence_quick_screen` | 120–151 | 32 |
| `refold_source_trace` | 152–445 | 294 |
| `raw_trajectory_trace` | 446–492 | 47 |

### 5.3 Candidate-Pool Maintenance

Low-confidence tiers are demoted, never deleted. If pH post-design finds a scaffold family more amenable to switching, related candidates can be pulled back from the trace tiers and put through structural filtering.

---

## 6. Scoring and Screening

### 6.1 Multi-Model Structural Filtering

Structural thresholds (BindCraft convention); all must be met for `structure_pass`:

| Metric | Threshold | What it measures |
|---|---|---|
| i_pTM | ≥ 0.70 | interface confidence |
| pTM | ≥ 0.55 | overall fold confidence |
| i_pAE | ≤ 0.35 | interface positional error |
| Unbound_Binder_pLDDT | ≥ 0.70 | fold confidence of the binder without its target |

The first three measure the complex's interface and fold confidence; `Unbound_Binder_pLDDT` measures the binder's fold confidence in the absence of the target. The two differ in input conditions.

### 6.2 Three-Species KD Predictor

Architecture: frozen protein-encoder embeddings → shared MLP → human / mouse / cyno heads, sharing the encoder and trunk while splitting the prediction head per species. The encoder is **ESM2-3B** (frozen).

**Ranking signals**:

- Primary: `cross_species_mean_prob` (mean probability across the three species)
- Secondary: `cross_species_min_prob` (weakest-species protection, guarding against collapse on any one species)

For P2, `cross_species_min_prob` acts as a gate rather than a ranking term alone: a candidate with a low mouse probability is demoted even if its mean is high.

**Training data distribution**:

| Species | Distribution | Source |
|---|---|---|
| Human | 923 usable records | Claude 90 + ProteinBase EGFR 833 |
| Mouse | binder 17 / non_binder 67 / not_tested 6 | Claude data only |
| Cyno | binder 10 / non_binder 74 / not_tested 6 | Claude data only |

The mouse and cyno heads are each trained on roughly 90 records with fewer than 20 positives, so their predictions are less reliable than the human head. P2 judgements must be cross-checked against conserved-epitope coverage (§3) rather than resting on this score alone.

### 6.3 pH Post-Design

**Principle**: the histidine side chain has pK<sub>a</sub> ≈ 6.0–6.5 — partially protonated and positively charged at pH 6.5, deprotonated and neutral at pH 7.4. Binding only becomes pH-dependent when an ionizable site titrates differently in the free and bound states, and the size of that dependence follows the proton-linkage relation

```
dG_site(pH)  = -RT · ln[ (1 + 10^(pKa_bound - pH)) / (1 + 10^(pKa_free - pH)) ]
ddG_switch   = dG(pH 6.5) - dG(pH 7.4)          negative = stronger binding at pH 6.5
```

**How many sites are needed**: one fully coupled site can shift binding by at most 2.303·RT·ΔpH = 1.23 kcal/mol across 6.5 → 7.4, about 8× in affinity. Losing detectable binding at pH 7.4 means roughly 100×, so designs carry **three or more coupled ionizable positions**, not one or two.

| Coupled sites | Ceiling | Affinity ratio |
|---|---|---|
| 1 | 1.23 kcal/mol | 8× |
| 2 | 2.46 kcal/mol | 63× |
| 3 | 3.69 kcal/mol | 502× |

**Titration on the target side**: PROPKA over 1YY9 chain A places EGFR **His409 inside the conserved epitope** with a free pK<sub>a</sub> of 6.24 — 35.5% protonated at pH 6.5 against 6.5% at pH 7.4. In the cetuximab complex that pK<sub>a</sub> falls to 4.96, so the antibody binding mode disfavours the protonated form and binds 0.2 kcal/mol *more weakly* at pH 6.5. Hotspot guidance copied from this epitope therefore biases towards anti-selectivity. A pH 6.5-ON design has to engage His409 in the opposite sense: present a carboxylate that stabilises the protonated state and raises pK<sub>a,bound</sub>. Glu472 (pK<sub>a</sub> 4.79) and Asp436 (pK<sub>a</sub> 3.82) lie 3.8 Å and 4.2 Å from the epitope and stay deprotonated across the window, so they serve as fixed negative partners for a binder-side histidine.

**Steps** (parents restricted to the de novo candidates of §5):

1. Select parents — all candidates in the `structure_pass_best` tier (rank 1–87)
2. Position selection — binder positions able to salt-bridge EGFR His409 take an Asp or Glu; binder histidines go adjacent to Glu472 or Asp436; anchor positions that hydrogen-bond or salt-bridge conserved epitope residues are excluded
3. Introduce three or more coupled ionizable positions per parent, enumerate combinations, deduplicate
4. Fold each variant, then run PROPKA on the complex, on the binder alone and on the target alone, and compute `ddg_switch` with the relation above
5. Gate on `ddg_switch ≤ -1.5 kcal/mol` (about 12×); survivors are re-scored through §6.1 and §6.2, and parent scores are not carried over
6. Composite ranking (each term converted to a percentile, then weighted):

```
final_score = 0.32 × predict          # affinity anchor
            + 0.38 × ddg_switch       # computed proton-linkage switch strength
            + 0.15 × ph_switch        # heuristic pH sensitivity, secondary
            + 0.10 × developability   # developability
            + 0.05 × mutation_burden  # penalty on edit count
```

The switch term outweighs the affinity term, matching P1 > P3; `predict` is retained as an anchor because added ionizable residues usually lower affinity, a cost the switch gain has to offset.

**Negative control**: running `scripts/ph_linkage.py` on 1YY9 with the Fab as binder returns `ddg_switch = -0.29 kcal/mol`, a 1.6× ratio, verdict FAIL — a known high-affinity, pH-insensitive binder is correctly rejected by the gate.

**Caveat**: PROPKA is an empirical predictor with a typical error near 0.5–1 pK<sub>a</sub> unit, and the numbers above come from a single crystal structure in a single conformation. APBS Poisson–Boltzmann or constant-pH molecular dynamics would tighten them at one to two orders of magnitude more compute.

---

## 7. Current Results and Candidate Pools

### 7.1 Leading De Novo Candidate

The most directly actionable tier is FreeBindCraft's `structure_pass_best`; the leading design is `egfr_bindcraft2_denovo_l97_7af9932070160527_candidate2` (selection_rank 1):

| Item | Value |
|---|---|
| Three-species probability | human 0.995 / mouse 0.754 / cyno 0.995, mean 0.915 |
| Structural metrics | i_pTM 0.90, pTM 0.92, i_pAE 0.12, Unbound_Binder_pLDDT 0.96 |

It meets both the structural and the three-species criteria. The RFdiffusion3 partial route has only completed a quick screen; refold and structural filtering come next.

### 7.2 Panel Selection

| Pool | Composition | Purpose |
|---|---|---|
| `global_top` | top N by composite score | scaffold-concentrated, probes the ceiling |
| `diverse_panel` | rebalanced across generation routes and backbone scaffolds | main body of the first assay panel, lowering whole-panel failure risk |

**First assay panel, suggested counts**: diverse_panel 60–80, global_top 20–40, high-affinity controls without pH design 8–12.

---

## 8. Submission

**Required CSV**:

| Column | Value |
|---|---|
| `name` | = `candidate_id` |
| `sequence` | = `seq` (10–250 aa) |
| `molecule_class` | `protein` |

**Ordering strategy**:

1. Highest pH-switch scores within `diverse_panel`
2. High scorers from `global_top`
3. Candidates with high `cross_species_min_prob`
