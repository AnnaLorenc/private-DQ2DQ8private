
📦 Python Implementation (PyMC v5)
import numpy as np
import pymc as pm
import arviz as az


def fit_public_model(
    Y_matrix: np.ndarray,      # shape (n_clones, n_individuals)
    depths: np.ndarray,        # shape (n_individuals,)
    min_shared: int = 3,
    draws: int = 2000,
    tune: int = 2000,
    target_accept: float = 0.9,
):
    """
    Bayesian occupancy model correcting for sequencing depth.

    Returns:
        posterior object
        probability each clone is shared in >= min_shared individuals
    """

    n_clones, n_individuals = Y_matrix.shape

    # Depth-based detection probability
    # small pseudo-frequency assumption
    f_bar = 1e-5
    p_detect = 1 - (1 - f_bar) ** depths
    p_detect = np.clip(p_detect, 1e-6, 1 - 1e-6)

    with pm.Model() as model:

        # Hyperprior for sharing probability
        a = pm.Exponential("a", 1.0)
        b = pm.Exponential("b", 1.0)

        # Clone-level true presence probability
        pi = pm.Beta("pi", alpha=a, beta=b, shape=n_clones)

        # Detection probability per individual
        p = pm.Data("p", p_detect)

        # Marginalized likelihood
        prob_obs = pi[:, None] * p[None, :]
        prob_obs = pm.math.clip(prob_obs, 1e-9, 1 - 1e-9)

        Y = pm.Bernoulli(
            "Y",
            p=prob_obs,
            observed=Y_matrix
        )

        trace = pm.sample(
            draws=draws,
            tune=tune,
            target_accept=target_accept,
            chains=4,
            cores=4,
        )

    # ---- Step 5: Compute probability shared >= N ----

    pi_samples = trace.posterior["pi"].values  # shape: chains, draws, clones
    pi_samples = pi_samples.reshape(-1, n_clones)

    prob_shared = []

    for j in range(n_clones):
        # For each posterior sample, simulate sharing count
        shared_counts = np.random.binomial(
            n_individuals,
            pi_samples[:, j]
        )
        prob = np.mean(shared_counts >= min_shared)
        prob_shared.append(prob)

    prob_shared = np.array(prob_shared)

    return trace, prob_shared

# 📊 Calling Public Clones (Step 5)
trace, prob_shared = fit_public_model(Y_matrix, depths, min_shared=5)

public_mask = prob_shared >= 0.95
public_clones = np.where(public_mask)[0]


# This yields:

# Posterior sharing probability per clone

# Robust classification

# 🧪 Step 6 — Validation & Sensitivity
# 6.1 Prior Sensitivity

# Re-run with:

a = pm.Exponential("a", 0.5)
b = pm.Exponential("b", 0.5)


# Compare posterior public clone set.

# 6.2 Compare with Rarefaction

# Perform:

# 200 subsampling runs

# Compare overlap with Bayesian calls

# Report:

# Jaccard index=
# ∣A∪B∣
# ∣A∩B∣
# 	​

# 6.3 Posterior Predictive Check
with model:
    ppc = pm.sample_posterior_predictive(trace)


# Verify:

# Observed sharing distribution

# Model-generated sharing distribution

# They should align.

# 🧬 What This Model Corrects

# ✔ Variable depth
# ✔ False negatives in shallow repertoires
# ✔ Overcalling due to deep repertoires
# ✔ Provides uncertainty

# ⚠ Important Assumption

# This version assumes:

# Detection probability depends only on depth

# Clone frequency distribution is similar across individuals

# We fit a hierarchical occupancy model correcting for depth-dependent detection probability. Public clonotypes were defined as those with posterior probability ≥0.95 of being present in at least N individuals.