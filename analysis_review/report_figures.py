"""Build report figures from verified cached tables, without rerunning analyses.

python analysis_review/report_figures.py
Numerical analysis and its independent checks remain in their original scripts.
"""
from pathlib import Path
import os
import tempfile

os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "venomics-mpl"))
os.environ["MPLBACKEND"] = "Agg"
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
DECISIONS = ROOT / "analysis_review/outputs/screening_decision_trends"
PANELS = ROOT / "analysis_review/outputs/shortlist_overlap"
FOCUSED = ROOT / "analysis_review/outputs/report_figures"
SPECIES = {"Gram-negative": "E. coli", "Gram-positive": "B. subtilis"}
COLORS = ["#225b78", "#bd6a34", "#287a65"]
CUTOFFS = [4, 8, 16, 32, 64, 128]


def style():
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 11,
                         "axes.spines.top": False, "axes.spines.right": False,
                         "pdf.fonttype": 42, "svg.fonttype": "none"})


def save_figure(fig, directory, name):
    directory.mkdir(parents=True, exist_ok=True)
    for ext in ["png", "pdf"]:
        fig.savefig(directory / (name + "." + ext), dpi=185, bbox_inches="tight")
    plt.close(fig)


def cutoff_axis(ax):
    ax.set(xscale="log", xlabel="Illustrative MIC cutoff (µM)")
    ax.set_xticks(CUTOFFS, labels=CUTOFFS)


def draw_cutoffs(ax, data):
    for (model, block), color in zip(data.groupby("model", sort=False), COLORS):
        ax.plot(block.cutoff_uM, block.pass_fraction * 100, "o-", label=model, color=color)
    cutoff_axis(ax)
    ax.set(ylabel="Peptides passing the rule (%)", ylim=(0, 105), title="Same cutoff, different candidate counts")
    ax.legend(frameon=False, fontsize=9)


def draw_matched(ax, data):
    for (target, species), color in zip(SPECIES.items(), COLORS):
        block = data[data.target.eq(target)]
        ax.plot(block.cutoff_uM, block.decision_disagreement_fraction * 100, "o-",
                label=f"{target} head ({species})", color=color)
    cutoff_axis(ax)
    ax.set(ylabel="CNN/LSTM pass–fail disagreement (%)", ylim=(0, 26), title="Same selected head, different decisions")
    ax.legend(frameon=False, fontsize=8)


def draw_votes(ax, row):
    values = np.array([row.all_positive_n, row.mixed_n, row.all_negative_n]) / row.n * 100
    ax.bar(["All positive", "Mixed", "All negative"], values, color=[COLORS[2], COLORS[1], "#a5b8c3"])
    for j, value in enumerate(values):
        ax.text(j, value + 1, f"{value:.1f}%", ha="center")
    ax.set(ylabel="Classifier-complete peptides (%)", ylim=(0, 75), title="Four classifiers often split their decisions")


def draw_mixed_votes(ax):
    patterns = pd.read_csv(DECISIONS / "classifier_vote_patterns.csv", dtype={"pattern": str})
    patterns = patterns[patterns.population.eq("all_classifier_complete")]
    counts = patterns.groupby(patterns.pattern.str.zfill(4).map(lambda p: p.count("1")))["count"].sum().loc[1:3]
    values = counts / counts.sum() * 100
    ax.bar(["1 positive", "2 positive", "3 positive"], values, color=["#dba37e", COLORS[1], "#8e4826"])
    for j, (count, value) in enumerate(zip(counts, values)):
        ax.text(j, value + .6, f"{value:.1f}%\n(n = {count:,})", ha="center", fontsize=9)
    ax.set(ylabel="Mixed-decision peptides (%)", ylim=(0, 47), title="How the mixed decisions split")


def decision_tables():
    counts = pd.read_csv(DECISIONS / "MIC_cutoff_counts.csv").query("population == 'all_common'")
    matched = pd.read_csv(DECISIONS / "matched_target_cutoff_decisions.csv").query(
        "comparison == 'DeepAMP_same_Gram_group' and population == 'all_paired'")
    votes = pd.read_csv(DECISIONS / "classifier_decision_counts.csv").query(
        "population == 'all_classifier_complete' and score_threshold == 0.5").iloc[0]
    return counts, matched, votes


def build_decision_panels():
    style()
    for draw, data, name in zip([draw_cutoffs, draw_matched], decision_tables()[:2],
                                ["concentration_cutoffs", "matched_target_decisions"]):
        fig, ax = plt.subplots(figsize=(7.5, 4.5), layout="constrained")
        draw(ax, data)
        save_figure(fig, PANELS, name)
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.7), layout="constrained")
    draw_votes(axes[0], decision_tables()[2])
    draw_mixed_votes(axes[1])
    fig.savefig(PANELS / "classifier_votes.svg", bbox_inches="tight")
    save_figure(fig, PANELS, "classifier_votes")


def build_classifier_sequence_figures():
    """Read the classifier-sequence analysis tables; keep inference separate."""
    out = ROOT / "analysis_review/outputs/classifier_sequences"
    if not (out / "sequence_profiles.csv").exists():
        return
    style()
    profiles = pd.read_csv(out / "sequence_profiles.csv")
    profiles = profiles[profiles.population.eq("all_complete") & profiles.grouping.eq("votes")].copy()
    profiles["votes"] = profiles["group"].astype(int)
    specs = [("length", "median", 1, "Median length (residues)"),
             ("KR_fraction", "mean", 100, "Mean K/R content (%)"),
             ("DE_fraction", "mean", 100, "Mean D/E content (%)"),
             ("hydrophobic_fraction", "mean", 100, "Mean hydrophobic content (%)"),
             ("contains_C", "mean", 100, "Sequences containing cysteine (%)"),
             ("entropy_bits", "mean", 1, "Mean residue diversity (bits)")]
    fig, axes = plt.subplots(2, 3, figsize=(12, 7), layout="constrained")
    vote_colors = ["#a5b8c3", "#dba37e", COLORS[1], "#8e4826", COLORS[2]]
    for ax, (feature, metric, scale, title) in zip(axes.flat, specs):
        block = profiles[profiles.feature.eq(feature)].sort_values("votes")
        ax.bar(block.votes, block[metric] * scale, color=vote_colors)
        ax.set(title=title, xticks=range(5), xlabel="Positive votes out of four")
    save_figure(fig, out, "sequence_profiles_by_votes")

    data = pd.read_csv(out / "source_length_matched_contrasts.csv")
    features = ["KR_fraction", "DE_fraction", "hydrophobic_fraction", "contains_C"]
    fig, axes = plt.subplots(1, 3, figsize=(13, 4.5), layout="constrained", sharey=True)
    comparisons = [("All negative", "Mixed"), ("All negative", "All positive"), ("Mixed", "All positive")]
    for ax, (left, right) in zip(axes, comparisons):
        block = data[data.left.eq(left) & data.right.eq(right)]
        for population, color, label in [("all_complete", COLORS[0], "Source + exact length"),
                                          ("mapped_equal_source_component", COLORS[1], "Also reduce group weight")]:
            values = block[block.population.eq(population)].set_index("feature").loc[features, "right_minus_left"] * 100
            ax.scatter(values, range(4), color=color, label=label, s=45)
        ax.axvline(0, color="#8999a1", lw=1)
        ax.set(title=f"{right} minus\n{left}", yticks=range(4),
               yticklabels=["K/R content", "D/E content", "Hydrophobic content", "Cysteine presence"],
               xlabel="Adjusted difference (percentage points)", ylim=(3.6, -.6), xlim=(-10, 16))
    axes[-1].legend(frameon=False, fontsize=8, loc="lower right")
    save_figure(fig, out, "matched_sequence_contrasts")

    transitions = pd.read_csv(out / "neighbor_group_transitions.csv")
    groups = ["All negative", "Mixed", "All positive"]
    matrix = transitions.pivot(index="anchor_group", columns="neighbor_group", values="fraction").loc[groups, groups] * 100
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.7), layout="constrained")
    axes[0].imshow(matrix, cmap="Blues", vmin=0, vmax=100)
    for i in range(3):
        for j in range(3):
            value = matrix.iloc[i, j]
            axes[0].text(j, i, f"{value:.1f}%", ha="center", va="center", color="white" if value > 60 else "#20313c")
    axes[0].set(xticks=range(3), yticks=range(3), xticklabels=groups, yticklabels=groups,
                xlabel="Neighbor's group", ylabel="Anchor's group", title="One-residue neighbors often keep their group")
    same = transitions[transitions.anchor_group.eq(transitions.neighbor_group)].set_index("anchor_group").loc[groups]
    axes[1].bar(np.arange(3) - .18, same.fraction * 100, width=.36, label="Observed neighbors", color=COLORS[0])
    axes[1].bar(np.arange(3) + .18, same.expected_same_group_source_length * 100, width=.36,
                label="Same source/length reference", color="#a5b8c3")
    axes[1].set(xticks=range(3), xticklabels=groups, ylabel="Same-group directed pairs (%)", ylim=(0, 100),
                title="Similarity exceeds the conditional reference")
    axes[1].legend(frameon=False, fontsize=8)
    save_figure(fig, out, "sequence_neighbor_groups")

    residues = pd.read_csv(out / "residue_composition.csv")
    matrix = residues.pivot(index="group", columns="amino_acid", values="mean_per_sequence_fraction").loc[groups] * 100
    fig, ax = plt.subplots(figsize=(12, 3.1), layout="constrained")
    plot = ax.imshow(matrix, cmap="YlGnBu", aspect="auto")
    ax.set(xticks=range(len(matrix.columns)), xticklabels=matrix.columns, yticks=range(3), yticklabels=groups,
           xlabel="Amino acid", title="Residue profiles overlap, with different average proportions")
    fig.colorbar(plot, ax=ax, label="Mean per-sequence residue content (%)")
    save_figure(fig, out, "residue_composition")


def paired_response_figure(data, group_column, groups, count_column, label):
    """Compare cached CNN/LSTM pair-response rates, with counts visible in labels."""
    fig, axes = plt.subplots(1, 2, figsize=(11.8, 4.5), layout="constrained", sharex=True)
    for ax, (target, species) in zip(axes, SPECIES.items()):
        labels = []
        for j, group in enumerate(groups):
            block = data[data.target.eq(target) & data[group_column].eq(group)].set_index("model")
            cnn, lstm = [100 * block.at[model, "twofold_change_fraction"] for model in ["CNN", "LSTM"]]
            ax.plot([cnn, lstm], [j, j], color="#a5b4be", lw=2)
            for model, value, color in [("CNN", cnn, COLORS[0]), ("LSTM", lstm, COLORS[1])]:
                ax.scatter(value, j, s=60, color=color, label=model if j == 0 else None, zorder=3)
            labels.append(label(group, int(block.iloc[0][count_column])))
        ax.set(title=f"{target} head\n(paper training: {species})", xlim=(0, 65),
               xlabel="Pairs with ≥2-fold MIC response to order (%)", yticks=range(len(groups)),
               yticklabels=labels, ylim=(len(groups) - .5, -.6))
        ax.legend(loc="lower left", frameon=False)
    return fig


def build_decision_figures():
    style()
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.7), layout="constrained")
    for ax, draw, data in zip(axes, [draw_cutoffs, draw_matched, draw_votes], decision_tables()):
        draw(ax, data)
        ax.title.set_fontsize(11)
    save_figure(fig, DECISIONS, "decision_trends")
    groups = {"all_available_pairs": "All composition pairs", "cysteine_free_length8_to39": "No C, length 8–39",
              "greedy_no_reused_component": "No reused group"}
    data = pd.read_csv(DECISIONS / "DeepAMP_composition_controlled_order.csv")
    fig = paired_response_figure(data, "population", groups, "n_pairs", lambda group, n: f"{groups[group]} ({n:,})")
    save_figure(fig, DECISIONS, "composition_controlled_order")


def build_focused_figures():
    style()
    strain = pd.read_csv(ROOT / "analysis_outputs/mic_variability_strain_labels.csv")
    labels = []
    for a, b, n in zip(strain.strain_a, strain.strain_b, strain.n):
        species = " ".join(a.split()[:2])
        pair = f"{species}: {a.removeprefix(species + ' ')} / {b.removeprefix(species + ' ')}"
        labels.append(f"{pair}\n(n = {n:,} paired numeric labels)")
    fig, ax = plt.subplots(figsize=(9.5, 5), layout="constrained")
    values = strain.fraction_over_fourfold * 100
    ax.barh(range(len(strain)), values, color=COLORS[0])
    for j, value in enumerate(values):
        ax.text(value + .35, j, f"{value:.1f}%", va="center")
    ax.set(yticks=range(len(strain)), yticklabels=labels, ylim=(len(strain) - .5, -.5), xlim=(0, 22),
           xlabel="Paired labels differing by more than fourfold (%)",
           title="Recorded MIC contrasts for the same sequence")
    save_figure(fig, FOCUSED, "strain_label_contrasts")

    data = pd.read_csv(ROOT / "analysis_review/outputs/model_scope_review/DeepAMP_order_low_prediction_sensitivity.csv")
    data = data.rename(columns={"export_target_label": "target"})
    group_column = "maximum_prediction_uM_each_sequence_and_model"
    groups = sorted(data[group_column].unique())
    fig = paired_response_figure(data, group_column, groups, "selected_pair_n",
                                lambda cap, n: f"All predictions ≤{cap:,.0f} µM\n({n:,} pairs)")
    save_figure(fig, FOCUSED, "order_response_low_prediction")

    overlap = pd.read_csv(ROOT / "analysis_outputs/experimental_training_overlap.csv")
    counts = overlap["Group (APEX training)"].value_counts().reindex(["CV", "Ind"], fill_value=0)
    fig, ax = plt.subplots(figsize=(7.1, 3.1), layout="constrained")
    start = 0
    for group, color, explanation in [("CV", COLORS[0], "CV corpus"), ("Ind", COLORS[1], "Ind corpus")]:
        n = int(counts[group])
        ax.barh([0], [n], left=start, height=.45, color=color, label=explanation)
        ax.text(start + n / 2, 0, f"{n} {group}", ha="center", va="center", color="white", weight="bold")
        start += n
    ax.set(xlim=(0, len(overlap)), xticks=range(0, len(overlap) + 1), yticks=[], ylim=(-.7, .7),
           xlabel="Numeric E. coli assay matches already present in APEX development/evaluation data",
           title=f"All {len(overlap)} numeric matches overlap the APEX corpus")
    ax.spines["left"].set_visible(False)
    save_figure(fig, FOCUSED, "experimental_label_exposure")


def main():
    build_decision_figures()
    build_decision_panels()
    build_focused_figures()
    build_classifier_sequence_figures()
    print("Rebuilt report figures from cached tables; no numerical tables changed.")


if __name__ == "__main__":
    main()
