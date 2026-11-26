# Linear-only / quadratic-only fits (morphing templates)

Source: Giacomo, whiteboard (photo `~/Desktop/DrellYan/65D327E4-74F1-454E-AAEA-F93E309AAA87.jpg`, sent on
Mattermost 2026-11-25), one operator. Generalisation to N operators (with w11) below: mine, checked numerically
(random SM, L, Q, M, k; residual 0 for full, lin-only, quad-only, 1 and 2 operators).
Status (2026-11-25): no combine model for this yet ("non c'è ancora una versione ma la faremo"); Giacomo will look at
it soon, so possibly just pick up his version. Low priority (exams, part time).

## Templates
Morphing templates = full predictions at a few Wilson-coefficient points, everything else = 0:
- w(0) = SM
- w(1)_i = SM + L_i + Q_i            (k_i = +1)
- w(-1)_i = SM - L_i + Q_i           (k_i = -1)
- w11_ij = SM + L_i + L_j + Q_i + Q_j + M_ij   (k_i = k_j = 1)

so L_i = (w1_i - wm1_i)/2, Q_i = (w1_i + wm1_i)/2 - SM, M_ij = w11_ij - w1_i - w1_j + SM.

## One operator (whiteboard)
Full (what the current model does):

    N(L+Q) = (1 - k^2) w(0) + k/2 (k+1) w(1) + k/2 (k-1) w(-1) = SM + k L + k^2 Q

Linear only: remove k^2 Q = k^2/2 (w(1) + w(-1) - 2 w(0)):

    N(L) = N(L+Q) - k^2 Q = w(0) + k/2 w(1) - k/2 w(-1) = SM + k L

Quadratic only (no interference), same logic:

    N(Q) = (1 - k^2) w(0) + k^2/2 w(1) + k^2/2 w(-1) = SM + k^2 Q

## N operators (with cross terms)
Coefficient of each template; sums over pairs i<j, "j != i" = all other operators:

| template | full: SM + sum k L + sum k^2 Q + sum k_i k_j M | lin only: SM + sum k L | quad only: SM + sum k^2 Q + sum k_i k_j M |
|---|---|---|---|
| w(0) | 1 - sum_i k_i^2 + sum_{i<j} k_i k_j | 1 | 1 - sum_i k_i^2 + sum_{i<j} k_i k_j |
| w1_i | k_i (k_i + 1)/2 - sum_{j!=i} k_i k_j | k_i/2 | k_i^2/2 - sum_{j!=i} k_i k_j |
| wm1_i | k_i (k_i - 1)/2 | -k_i/2 | k_i^2/2 |
| w11_ij | k_i k_j | 0 | k_i k_j |

- Lin only needs no w11 (mixed terms are O(1/Lambda^4)), so it works on the 1D cards too.
- "Quad only" here keeps the mixed terms M_ij (they are dim-6 squared too); drop the w11 column and the k_i k_j
  terms for squares only.

## Caveats
- Lin only: SM + k L can go negative in a bin for large |k| (nothing like k^2 Q keeps it positive) -> combine fails
  on negative yields (see memory: EFT negative-yield boundary). Keep the scan ranges near the linear regime.
- The coefficient of wm1 is negative: fine for a sum of scaled processes, but the per-bin total must stay > 0.
- Physics use: lin only = O(1/Lambda^2), standard truncation; lin vs full shows the size of the dim-6^2 terms
  (stand-in for the missing dim-8). Quad only is a diagnostic, not a physical truncation.

## Implementation (to do)
- The morphing model (Giacomo's, `createWS.py` / physics model in tools) builds these coefficients for the full case;
  a lin-only option = different coefficient functions per template (table above), same cards.
- Check first whether Giacomo already has it (he said he will look at it).
