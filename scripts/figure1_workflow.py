"""Figure 1: study workflow (revised; no synthetic data). Usage: python figure1_workflow.py <figures_dir>"""
import sys, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

FIG = sys.argv[1]
INK, INK2, EDGE = '#0b0b0b', '#52514e', '#9a9994'
FILL = {'cohort': '#f0efec', 'data': '#cde2fb', 'desc': '#eaf3ec', 'model': '#fdeee6', 'ext': '#efedf9'}
plt.rcParams.update({'font.family': 'DejaVu Sans', 'pdf.fonttype': 42})
fig, ax = plt.subplots(figsize=(12, 8.6)); ax.set_xlim(0, 106); ax.set_ylim(0, 86); ax.axis('off')

def box(x, y, w, h, title, body, kind):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle='round,pad=0.4,rounding_size=1.2', fc=FILL[kind], ec=EDGE, lw=0.9))
    ax.text(x + w / 2, y + h - 2.2, title, ha='center', va='top', fontsize=10, fontweight='bold', color=INK)
    ax.text(x + w / 2, y + h - 6.2, body, ha='center', va='top', fontsize=8.6, color=INK2, linespacing=1.45)

def arrow(x1, y1, x2, y2):
    ax.annotate('', xy=(x2, y2), xytext=(x1, y1), arrowprops=dict(arrowstyle='-|>', color=INK2, lw=1.1, shrinkA=2, shrinkB=2))

box(22, 72, 62, 12.5, 'Prospective cohort (2019–2021)',
    '42 patients with an MRI-visible lesion and radical prostatectomy\nIndex lesion: MRI region of interest matched to whole-mount pathology\nEndpoint: BCR (PSA > 0.2 ng/mL, confirmed) vs no evidence of disease', 'cohort')
box(6, 52, 42, 13, 'Transcriptome of the index lesion',
    'Ion AmpliSeq Transcriptome (FFPE RNA)\n35 patients: 11 BCR, 24 NED', 'data')
box(58, 52, 42, 13, 'T2-weighted MRI radiomics',
    'PyRadiomics features from the index lesion\n34 patients: 11 BCR, 23 NED', 'data')
arrow(42, 72, 27, 65.4); arrow(64, 72, 79, 65.4)
box(6, 31, 42, 15.5, 'Descriptive analyses',
    'Differential expression (DESeq2, FDR)\nGene set enrichment (Hallmark,\nphenotype permutation)', 'desc')
box(58, 31, 42, 15.5, 'Descriptive analyses',
    'Between-group comparisons and\ncorrelation differences\n(Benjamini–Hochberg FDR)', 'desc')
arrow(27, 52, 27, 46.9); arrow(79, 52, 79, 46.9)
box(6, 8, 42, 17.5, 'Explanatory models (descriptive only)',
    'XGBoost fitted on all patients + SHAP\ntranscriptomic, radiomic, combined\nfeature stability across bootstrap\nrefits (Figure 4)', 'model')
box(58, 8, 42, 17.5, 'Predictive performance',
    'Repeated stratified 5-fold CV;\nall supervised steps inside training\nfolds; held-out real patients only\n(Supplementary Table S5)', 'model')
for x in (27, 79): ax.plot([x, x], [31, 28.6], color=INK2, lw=1.1)
ax.plot([27, 79], [28.6, 28.6], color=INK2, lw=1.1)
arrow(27, 28.6, 27, 26.2); arrow(79, 28.6, 79, 26.2)
for ext in ('pdf', 'png'):
    fig.savefig(f'{FIG}/Figure1.{ext}', dpi=300, bbox_inches='tight')
print('Figure 1 written')
