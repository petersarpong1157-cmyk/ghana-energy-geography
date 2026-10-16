# Ghana Energy Geography

Research website for **Renewable Energy Opportunity and the Geography of Infrastructure Disadvantage: Settlement-Level Evidence from Ghana**.

**Authors:** Peter Sarpong and Dora Konadu  
**Author details:** Peter Sarpong — Independent Researcher; Dora Konadu — graduate of the University of Education, Winneba, Ghana (2023). Her current institutional affiliation has not been specified.

## Website
The landing page is in [index.html](index.html), designed to run as a static GitHub Pages site. It includes responsive styling, descriptive data visualizations, research methods, key numerical results, and scientific limitations.

**Website address:** https://petersarpong1157-cmyk.github.io/ghana-energy-geography/ (GitHub Pages was enabled by the repository owner; external live availability has not yet been independently confirmed.)

**Reproducibility status:** See [REPRODUCIBILITY.md](REPRODUCIBILITY.md) for verified tests and remaining caveats. Use [CITATION.cff](CITATION.cff) for repository citation metadata.

## Empirical replication files

The original 11 R scripts and verified Phase 6–9 output tables are now available in [analysis/](analysis/README.md). The upstream source datasets are linked rather than redistributed publicly pending rights/provenance review. This is a **partially independently validated research archive**, not a certified end-to-end replication.

## Study summary
- Source records: 25,342 settlement records from the World Bank Ghana DRE Atlas.
- Main complete-case model: 21,582 observations.
- Principal findings: modeled relative wealth is negatively associated with distance to transmission lines and roads, conditional on model covariates and regional controls.
- Solar potential association is sensitive to geographic specification and spatial inference bandwidth.

## Interpretation safeguards
- **No causal interpretation.** Proximity to transmission lines is not verified electricity access.
- The Relative Wealth Index is modeled and could share upstream geospatial information with infrastructure measures.
- Residual spatial dependence and uncertainty about boundary provenance remain.
- Solar potential represents long-term modeled resource potential, not installed PV capacity or electricity generation.

## Data source
[World Bank Ghana Distributed Renewable Energy (DRE) Atlas](https://energydata.info/dataset/ghana-distributed-renewable-energy-dre).

## Research outputs
The final paper, data/code archive, and permanent SSRN/DOI links are **not** included here unless independently verified and explicitly approved for public release.

## Site publication
In repository **Settings → Pages**, set **Build and deployment → Deploy from a branch**, choose **main** and **/(root)**, then save. The intended address is https://petersarpong1157-cmyk.github.io/ghana-energy-geography/ .

## AI assistance
ChatGPT was used for research planning, coding support, interpretation, manuscript drafting and language editing. The researcher is responsible for verification and conclusions.
