"""Draws Figure 1 (pathway identification flow). The counts are those of the
pathway selection (01_pathway_selection.R; Supplementary Table S2b).

Hinz et al. Pathway-restricted polygenic scores for energy-substrate metabolism
and obesity (FoCus cohort). See README for inputs and run order.

Output: <TABLES_OUT_DIR> (default: ./tables_figures_output)
"""
import os
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
plt.rcParams["font.family"] = "DejaVu Sans"
fig, ax = plt.subplots(figsize=(10.2, 12.4))
ax.set_xlim(0, 100); ax.set_ylim(0, 124); ax.axis("off")
EDGE="#222222"; BLUE="#EEF3F7"; GREY="#F5F5F5"; RED="#B22222"
def box(x,y,w,h,text,fill=GREY,fs=10.5,bold_first=False):
    ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle="round,pad=0.25,rounding_size=0.8",fc=fill,ec=EDGE,lw=1.4))
    L=text.split("\n")
    if bold_first:
        ax.text(x+w/2,y+h-2.0,L[0],ha="center",va="top",fontsize=fs,weight="bold")
        ax.text(x+w/2,y+(h-3.0)/2,"\n".join(L[1:]),ha="center",va="center",fontsize=fs-0.5,linespacing=1.35)
    else:
        ax.text(x+w/2,y+h/2,text,ha="center",va="center",fontsize=fs,linespacing=1.35)
def arrow(x1,y1,x2,y2,dashed=False,col=EDGE):
    ax.add_patch(FancyArrowPatch((x1,y1),(x2,y2),arrowstyle="-|>",mutation_scale=16,lw=1.6,color=col,
                 linestyle=(0,(4,3)) if dashed else "-"))
def excl(x,y,text):
    ax.plot([x,x+8],[y,y],color=RED,lw=1.3); ax.text(x+8.5,y,text,ha="left",va="center",fontsize=9.5,color=RED)

ax.text(50,122,"Energy-substrate-metabolism pathway identification",ha="center",va="center",fontsize=14)
box(5,103,90,15,"Pre-specified selection criteria\n"
    "Include: macronutrient flux, storage, mobilization or oxidation; central nodes of energy\n"
    "metabolism; direct hormonal/metabolic regulation of substrate flux.\n"
    "Exclude: xenobiotic/secondary metabolism, non-energetic biosynthesis, pure signal\n"
    "transduction, disease-specific pathways.  Formal and size thresholds: Table S1",fill="#FFFFFF",bold_first=True)

# row B
box(5,84,38,14,"KEGG (Release 119.0)\nBRITE core categories: carbohydrate 15,\nenergy 8, lipid 17, amino acid 14\nn = 44",fill=BLUE,bold_first=True)
box(57,84,38,14,"Reactome (release 97)\nCandidates screened in 4\ncorresponding top-level branches\nn = 43",fill=BLUE,bold_first=True)
arrow(24,83.6,24,74.4); excl(24,79,"excluded: n = 18")
arrow(76,83.6,76,74.4); excl(76,79,"excluded: n = 18")
# row C
box(5,63,38,11,"Retained after criteria: n = 26\n+ gap search, independent of\nBRITE category: n = 9")
box(57,63,38,11,"Reactome candidates after criteria\nn = 25\nTwo-zone Jaccard model (J = 0.3)")
arrow(24,62.6,24,53.4)
arrow(70,62.6,70,53.4); arrow(87,62.6,87,53.4)
# row D
box(5,42,38,11,"35 KEGG pathways\n(6 of these enriched with\nnon-overlapping Reactome genes)")
box(57,42,24,11,"Zone 1: J < 0.3\nindependent\nn = 19")
box(83,42,12,11,"Zone 2:\nJ ≥ 0.3\nn = 6",fs=10)
ax.add_patch(FancyArrowPatch((86,41.6),(43.4,45),arrowstyle="-|>",mutation_scale=14,lw=1.3,color="#555555",
             linestyle=(0,(4,3)),connectionstyle="arc3,rad=-0.25"))
ax.text(90.5,36.6,"non-overlapping genes\nmerged into KEGG",ha="center",va="center",fontsize=8.8,color="#555555")
# row E
box(20,18,60,12,"54 pathway-level gene sets\n35 KEGG + 19 Reactome\ngene identifiers harmonized to HGNC symbols",fill=BLUE,bold_first=True)
arrow(24,41.6,38,30.4); arrow(69,41.6,60,30.4)
box(20,1,60,11,"SNP extraction, QC and LD pruning (r² 0.2)\n14–1 354 SNPs per pathway (median 189)\n→ pathway SNP-sum-scores (Figure 2)")
arrow(50,17.6,50,12.4)
OUT = Path(os.environ.get("TABLES_OUT_DIR", "tables_figures_output")); OUT.mkdir(parents=True, exist_ok=True)
out=str(OUT / "Figure1_pathway_selection")
plt.savefig(out+".png",dpi=300,bbox_inches="tight",facecolor="white")
plt.savefig(out+".pdf",bbox_inches="tight",facecolor="white")
plt.savefig(out+".tiff",dpi=300,bbox_inches="tight",facecolor="white")
