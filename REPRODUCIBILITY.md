# Reproducibility and scientific verification

**Project:** Renewable Energy Opportunity and the Geography of Infrastructure Disadvantage: Settlement-Level Evidence from Ghana  
**Researcher:** Peter Sarpong, Independent researcher  
**Status:** Working-paper research; the repository currently hosts a public research landing page, not a complete independent replication archive.

## Primary data
World Bank Ghana Distributed Renewable Energy (DRE) Atlas: https://energydata.info/dataset/ghana-distributed-renewable-energy-dre

- Settlement records in source: **25,342**
- Main complete-case regression observations: **21,582**
- Common polygon-harmonization analytical sample: **21,576**
- Source observations missing modeled relative wealth: **3,735**

Do not confuse proximity to transmission lines with verified electricity access or modeled solar potential with observed electricity production.

## Numerical verification
- Region fixed-effects principal regression was independently reproduced from the source CSV in Python.
- Reported Phase 6 and Phase 7 numerical output comparisons agree; Phase 9 common-sample harmonization results were reconciled.
- This repository now includes [all eleven archived R phase scripts and the available Phase 6–9 numerical tables](analysis/README.md). It does **not** yet include all upstream source inputs or a fully independent end-to-end rerun of the spatial HAC and polygon pipeline.

## Outstanding qualifications
1. The supplied district boundary archives contain ArcGIS processing timestamps from 2023; these do **not** establish legal boundary vintage or authoritative provenance.
2. The modeled Relative Wealth Index uses remotely sensed and connectivity-related information. Exact overlap with regressors remains unresolved.
3. Spatial autocorrelation remains after region and district fixed effects.
4. Missing wealth measurements may entail sample-selection bias; inverse-probability weighting addresses only modeled observables under stated assumptions.
5. Solar associations depend on fixed-effects choice and spatial correlation bandwidth.
6. A fully independent rerun of the Conley and polygon-harmonized pipelines remains outstanding.

## Planned archive organization
The scripts and selected result tables are published under `analysis/`. Do not upload third-party datasets or other materials until their licensing and sharing permissions are confirmed.

## Manuscript link
A verified SSRN abstract URL/DOI has **not yet been supplied**. When confirmed, add it here and update the website. Neither this page nor the GitHub repository implies SSRN acceptance or peer review.

## AI assistance
ChatGPT contributed to analysis planning, coding support, interpretation, manuscript drafting, and editing. The researcher remains responsible for verification and claims.
