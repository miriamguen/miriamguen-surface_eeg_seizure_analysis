"""To create a plot illustrating the Bayesian measures (Probability of Direction (PD), Region of Practical Equivalence (ROPE), and MAP-based p-value), we need to follow these steps:

Generate Synthetic Data: Create synthetic data to represent posterior distributions for illustration purposes.
Plot Each Measure: Create subplots for each measure, including Probability of Direction, MAP-based p-value, and ROPE.
Here's a Python script using Matplotlib and Seaborn to generate the required plots:

"""

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


# Function to create the Bayesian measures plot
def plot_bayesian_measures(save_path):
    # Generate synthetic data
    np.random.seed(0)
    posterior = np.random.normal(1, 0.5, 100000)
    prior = np.random.normal(0, 0.5, 100000)

    # Probability of Direction (pd)
    pd = (posterior > 0).mean() * 100

    # MAP-based p-value
    density_posterior = sns.kdeplot(posterior, bw_adjust=0.5).get_lines()[0].get_data()

    density_prior = sns.kdeplot(prior, bw_adjust=0.5).get_lines()[0].get_data()

    plt.close()
    map_posterior = density_posterior[0][np.argmax(density_posterior[1])]
    density_at_0_posterior = density_posterior[1][np.abs(density_posterior[0]).argmin()]
    map_p_value = density_at_0_posterior / np.max(density_posterior[1])

    # ROPE
    rope_min, rope_max = -0.1, 0.1
    rope_full = ((posterior > rope_min) & (posterior < rope_max)).mean() * 100
    hdi_95 = np.percentile(posterior, [2.5, 97.5])
    rope_95 = (
        (posterior > rope_min)
        & (posterior < rope_max)
        & (posterior > hdi_95[0])
        & (posterior < hdi_95[1])
    ).mean() * 100

    # Create subplots
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))

    # Plot PD
    sns.histplot(
        posterior, bins=100, kde=True, ax=axes[0], color="royalblue", stat="density"
    )
    for patch in axes[0].patches:
        if patch.get_x() < 0:
            patch.set_facecolor("orange")
        else:
            patch.set_facecolor("royalblue")
        patch.set_alpha(0.5)

    axes[0].axvline(0, color="orange", linestyle="--")
    # axes[0].fill_betweenx(y=[0, axes[0].get_ylim()[1]], x1=-4, x2=0, color='red', alpha=0.3)
    axes[0].set_title(f"Probability of Direction (pd): {pd:.1f}%")
    axes[0].set_xlabel("Value")
    axes[0].set_ylabel("Density")

    # Plot MAP-based p-value
    sns.histplot(
        posterior, bins=100, kde=True, ax=axes[1], color="royalblue", stat="density"
    )
    axes[1].axvline(map_posterior, color="k", linestyle="--", label="Empirical MAP")
    axes[1].axvline(0, color="orange", linestyle="--", label="Null Hypothesis")
    axes[1].legend()
    axes[1].set_title(f"MAP-based p-value: {map_p_value:.3f}")
    axes[1].set_xlabel("Value")

    # Plot ROPE
    sns.histplot(
        posterior, bins=100, kde=True, ax=axes[2], color="royalblue", stat="density"
    )
    # axes[2].axvspan(rope_min, rope_max, color='red', alpha=0.3, label='ROPE')
    for patch in axes[2].patches:
        if (patch.get_x() > rope_min) and (patch.get_x() < rope_max):
            patch.set_facecolor("orange")
        else:
            patch.set_facecolor("royalblue")
        patch.set_alpha(0.5)

    axes[2].axvline(0, color="orange", linestyle="--", label="Null Hypothesis")
    axes[2].legend()
    axes[2].set_title(f"ROPE (Full): {rope_full:.1f}%, ROPE (95%): {rope_95:.1f}%")
    axes[2].set_xlabel("Value")

    plt.tight_layout()
    plt.savefig(save_path, format="svg", dpi=300)
    plt.close()


# Example usage
plot_bayesian_measures("bayesian_measures.svg")
