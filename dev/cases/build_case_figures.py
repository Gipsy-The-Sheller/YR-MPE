"""Assemble the three case-study figures from the rendered UI assets."""

from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

plt.rcParams["font.family"] = "Arial"

HERE = Path(__file__).resolve().parent
ASSETS = HERE / "assets"
OUT = HERE / "figures"
OUT.mkdir(parents=True, exist_ok=True)

CASES = {
    "case_distance.png": [
        ("A", "ITS2 ML distance matrix", "F_its2_distance_matrix.png",
         "Eight-taxon subset; values are ML distances under TIM2+F+G4."),
    ],
    "case_species_tree.png": [
        ("A", "ASTRAL-IV species tree", "A_astral_tree.png",
         "Species tree from the 150-gene-tree ILS data set."),
        ("B", "gCF/qCF simplex", "ils_simplex.png",
         "Each point is an internal split; colour is gDFP."),
        ("C", "TreeRecs reconciliation", "C_treerecs.png",
         "Gene/species reconciliation with duplication and loss symbols."),
    ],
    "case_selection.png": [
        ("A", "Codon Z-test table", "A_selection_table.png",
         "Pairwise NG86 estimates and selection calls."),
        ("B", "MEGA-style matrix", "B_selection_matrix.png",
         "Upper triangle: two-sided p; lower triangle: dN/dS."),
    ],
    "manual_msc_simplex.png": [
        ("A", "gCF and ASTRAL-IV qCF simplex", "ils_simplex.png",
         "Each point is one internal split; gDFP is encoded by colour."),
        ("B", "Bipartition membership matrix", "ils_bipartition.png",
         "Taxa as rows, splits as columns; blue marks the smaller side."),
    ],
    "case_seqmatrix_colors.png": [
        ("A", "NCBI download: blue", "D_seqmatrix_ncbi_blue.png",
         "Cells filled from accession downloads are blue."),
        ("B", "Duplicate accession: red", "E_seqmatrix_duplicate_red.png",
         "The same accession appears in two cells; both are red."),
    ],
    "case_single_gene.png": [
        ("A", "ITS2 alignment", "A_its2_alignment.png",
         "25 sequences after MAFFT with direction correction and TrimAl."),
        ("B", "ModelFinder", "B_its2_modelfinder.png",
         "BIC-ranked substitution models; TIM2+F+G4 selected."),
        ("C", "IQ-TREE 3 tree", "C_its2_tree.png",
         "Maximum-likelihood tree with UFBoot support, shown in IcyTree."),
    ],
    "case_sanger.png": [
        ("A", "SeqMatrix: 12 taxa × 3 amplicons", "A_seqmatrix_sanger.png",
         "One missing ITS amplicon is left blank."),
        ("B", "MAFFT strand correction", "B_mafft_adjustdirection.png",
         "Advanced > Direction adjustment > --adjustdirection."),
        ("C", "Concatenated ML tree", "C_sanger_tree.png",
         "IQ-TREE 3 tree shown in IcyTree."),
    ],
    "case_organelle.png": [
        ("A", "SeqDBG: comparative chloroplast genomics", "A_seqdbg_organelle.png",
         "Gene statistics and adjacency Sankey for four annotated plastomes."),
        ("B", "Organelle species tree", "B_organelle_tree.png",
         "Tree from the linked organelle CDS workflow, shown in IcyTree."),
    ],
    "case_infcat_setup.png": [
        ("A", "PhyloBayes-MPI: infinite CAT", "A_phylobayes_infcat.png",
         "GTR + discrete gamma + infinite CAT; two independent chains."),
    ],
    "case_mrbayes.png": [
        ("A", "MiniTracer: MrBayes trace and density", "B_minitracer_trace_density.png",
         "Post-burn-in trace, marginal density, ESS and 95% HPD."),
        ("B", "MiniTracer: MrBayes two-chain comparison", "C_minitracer_compare_psrf.png",
         "Chain 1 vs chain 2 with the Gelman-Rubin PSRF (R-hat)."),
    ],
    "case_bayes.png": [
        ("A", "PhyloBayes-MPI: infinite CAT", "A_phylobayes_infcat.png",
         "GTR + discrete gamma + infinite CAT; two independent chains."),
        ("B", "MiniTracer: trace and density", "B_minitracer_trace_density.png",
         "Post-burn-in trace, marginal density, ESS and 95% HPD."),
        ("C", "MiniTracer: two-chain comparison", "C_minitracer_compare_psrf.png",
         "Chain 1 vs chain 2 with the Gelman–Rubin PSRF (R-hat)."),
    ],
}


def draw_case(filename: str, panels) -> None:
    count = len(panels)
    fig, axes = plt.subplots(1, count, figsize=(4.9 * count, 5.7),
                             facecolor="white")
    if count == 1:
        axes = [axes]
    for axis, (letter, title, asset, _note) in zip(axes, panels):
        image = mpimg.imread(ASSETS / asset)
        axis.imshow(image, interpolation="lanczos")
        axis.set_axis_off()
        axis.set_anchor("N")
        axis.add_patch(Rectangle((0, 0), 1, 1, transform=axis.transAxes,
                                 fill=False, edgecolor="#d1d5db",
                                 linewidth=0.8, clip_on=False))
        axis.set_title("%s  %s" % (letter, title), loc="left", pad=10,
                       fontsize=12, fontweight="bold", color="#111827")
    fig.subplots_adjust(left=0.012, right=0.988, top=0.93, bottom=0.03,
                        wspace=0.06)
    fig.savefig(OUT / filename, dpi=220, bbox_inches="tight", pad_inches=0.04)
    fig.savefig(OUT / filename.replace(".png", ".pdf"), bbox_inches="tight",
                pad_inches=0.04)
    plt.close(fig)
    print(OUT / filename)


def main() -> None:
    for filename, panels in CASES.items():
        draw_case(filename, panels)


if __name__ == "__main__":
    main()
