# Figure Generation Metadata

## Figure 1: Study Workflow Diagram

**Tool Used:** Google Gemini (Image Generation) / Nanobanana
**Date:** 2024-12-06
**Purpose:** To visualize the study design and workflow for the manuscript.

Note: the date above records when the prompt text was drafted. The `plots/publication/figure1.jpeg` file included in this package may be regenerated later using the same prompt.

### Prompt Used:
```
Please generate a publication-quality research flowchart for an SCI journal article.

IMPORTANT: 
- DO NOT include any title or figure number at the top of the image (the caption will be added separately below the figure)
- Use LANDSCAPE orientation with 4:3 aspect ratio (width > height)
- Arrange modules in a more horizontal layout to fit the landscape format

【Flowchart Structure】(5 main modules, arranged to fit 4:3 landscape)

=== Module 1: Data Sources (Top row, spread horizontally) ===
Left side:
- GEO Database: GSE14520 (HCC cohort, n=445, Tumor: 225, Non-tumor: 220)
- TCGA Database: TCGA-LIHC (External validation, n=369)

Right side:
- Ferroptosis Gene Database: FerrDb V2 (108 genes)
- TCM Database: TCMSP - Artemisia capillaris targets (42 targets)
  - 6 active compounds: Scoparone, Capillarisin, Chlorogenic acid, Quercetin, Kaempferol, Isorhamnetin

=== Module 2: DEG & Intersection Analysis ===
- GSE14520: Tumor vs Non-tumor → DEGs (1,391 genes)
- Three-set Venn diagram: DEGs ∩ Ferroptosis genes ∩ TCM targets
- → 3 Hub genes: ACSL4, TFRC, NQO1
- GO/KEGG enrichment analysis

=== Module 3: Prognostic Model Construction ===
- Univariate Cox regression (p<0.1): 24 genes → LASSO-Cox → 15-gene Risk Score model
- Performance: C-index=0.717, 1-yr AUC=0.742, 3-yr AUC=0.765, Log-rank p=3.1×10⁻⁸
- Clinical utility: Nomogram + DCA + Forest plot

=== Module 4: Functional & Mechanism Analysis (can be side-by-side) ===
Left panel:
- Immune infiltration (ssGSEA, 25 cell types): MDSC↑(r=0.48), M2 Macrophage↑(r=0.42)
- Immune checkpoint analysis (14 genes)
- Drug sensitivity (GDSC): 5-FU (p=0.001), Sorafenib (p=0.06)

Right panel:
- Network pharmacology: Herb-Compound-Target-Disease network

=== Module 5: Molecular Docking Validation (Dual Methods, bottom row) ===
Traditional Docking (left):
- CB-Dock2 (AutoDock Vina), 18 combinations (3 proteins × 6 ligands)
- All effective binding (Vina < -5.0 kcal/mol)
- Best: NQO1 + Isorhamnetin (-9.1 kcal/mol)

AI-based Validation (right):
- NVIDIA Boltz-2
- Best: ACSL4 + Scoparone (ipTM=0.833)
- → Convergent results support therapeutic potential

【Design Requirements】
1. Layout: LANDSCAPE 4:3 aspect ratio (e.g., 2400×1800 pixels or 170mm×127mm)
2. NO title/figure number at top - leave clean for journal caption
3. Arrange modules more horizontally to utilize landscape format
4. Color scheme: 
   - Data sources: Light Blue (#E3F2FD)
   - Analysis steps: Light Green (#E8F5E9)
   - Validation: Light Orange (#FFF3E0)
   - Key results: Light Red (#FFEBEE)
5. Icons: Add small icons for databases, methods, key findings
6. Highlight hub genes (ACSL4, TFRC, NQO1) in bold
7. Style: Clean, professional, SCI journal quality
8. Resolution: 300 dpi minimum
9. Background: White, no watermarks
10. Clear arrows connecting modules in logical flow
```

### Note on AI Usage:
This figure was generated using AI assistance to visualize the study design based on the specific prompt above. The scientific logic, data analysis, and workflow steps were conceived and designed entirely by the authors. All other figures (Figures 2-6) were generated using R scripts provided in the `scripts/` directory of this repository.
