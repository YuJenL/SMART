# SMART
SMART written in R and R GUI has been developed as user-friendly software for integrated analysis of metabolomics data. 

SMART streamlines the complete analysis flow from initial data preprocessing to advanced downstream data analysis, consisting of the following eight modules. 
1.	Data Import: Analyze different data file formats (e.g., .raw, .d, and mzXML).
2.	Data Visualization: Visually represent various types of data features (e.g., total ion chromatogram (TIC) and mass spectra).
3.	Peak Analysis and annotation: Implement peak analysis for both untargeted and targeted data and peak annotation.
4.	Data Preprocessing: Perform data transformation and standardization.
5.	Quality Control: Conduct quality control for samples and peaks.
6.	Batch Effect Analysis: Explore batch effects (e.g., known experimental conditions, unknown latent groups (LGs), or hidden substructures).
7.	Statistical Analysis: Perform association analysis (ANCOVA), classification analysis (PLS/PLS-DA), and integrative omics pathway analysis (IOPA; academic use only, see [License and third-party components](#license-and-third-party-components)).
8.	Post-analysis: Execute peak identification and concentration calibration.

## Notice (November 2026): KEGG-derived data files removed

The KEGG-derived data files used by the IOPA module (`hsaSPIA_reversible(cumulative_sum).RData`
and `KEGG_pathway_class.Rdata`) are no longer distributed in this repository, and they have been
removed from the repository history. KEGG data are copyrighted by Kanehisa Laboratories, and their
license terms do not permit redistribution without a license.

- The IOPA method is unchanged and is described in the SMART 2.0 paper (see Citations below).
- All other SMART modules are unaffected.
- Because the repository history was rewritten, commit IDs created before this change (including the `v2.0` tag) have changed.
- Academic users who wish to run IOPA may contact the corresponding author (see Correspondence below).

## License and third-party components

SMART's own source code is released under the GNU General Public License v3.0 (see `LICENSE`).
The GPL-3.0 applies only to code written by the SMART team. The following third-party
components and data used by the IOPA module are **not** covered by the GPL-3.0 and remain
subject to their own license terms:

- **SPIA** (Tarca AL, et al. A novel signaling pathway impact analysis. *Bioinformatics* 2009;25(1):75–82).
  `spia_mod()` in `4_IOPA.R` is adapted from the SPIA Bioconductor package, and the IOPA module
  also calls SPIA functions. SPIA is freely available to academic researchers for research and
  teaching purposes only; any other use requires a license from Advaita Bioinformatics.
- **KEGG pathway data**: the IOPA module uses pathway data derived from the KEGG database
  (© Kanehisa Laboratories): `hsaSPIA_reversible(cumulative_sum).RData` and `KEGG_pathway_class.Rdata`.
  Because of KEGG license terms, **these data files are not included in this distribution**.
  Without them, all other SMART modules work normally, and the IOPA module shows a notice instead
  of running. KEGG data are provided for academic use only; non-academic use requires a license
  from Pathway Solutions (https://www.pathway.jp).

**Commercial or other non-academic users** must obtain the appropriate licenses from
Advaita Bioinformatics and Pathway Solutions before using the IOPA module.
All other SMART modules do not depend on these components.

## ☁️ iMAS Cloud Platform

SMART 2.0 is freely available for academic use (SMART's own code is released under GPL-3.0;
see [License and third-party components](#license-and-third-party-components)).
However, setting up the required computational environment (e.g., Linux/WSL, Docker, and
high-memory instances) to process large-scale untargeted or targeted MS datasets can be
resource-intensive.

The **iMAS (Intelligent Metabolomics Analytical System) Team** offers a managed, ready-to-use
cloud platform hosted on a national high-performance computing facility in Taiwan, available
in two editions: **SMART 1.1 (Basic)** and **SMART 2.0 (Professional)**.

* **No IT setup:** no installation or container configuration required.
* **High performance:** computing resources for large MS datasets, as specified in your
  iMAS subscription agreement.
* **Support:** maintenance and technical assistance by the iMAS Team, as specified in your
  iMAS subscription agreement.

Service scope, computing resources, and service levels are governed by the iMAS subscription
agreement. In case of any difference, the agreement prevails.

> **Note:** The IOPA (pathway analysis) module is **not** available on the cloud platform
> because of third-party license terms (see Notice above).

👉 For subscription plans, please contact: lyj626@stat.sinica.edu.tw

*Citations:                                                                                              
Yu-Jen Liang, Yu-Ting Lin, Chia-Wei Chen, Chien-Wei Lin, Kun-Mao Chao, Wen-Harn Pan and Hsin-Chou Yang (2016/05). SMART: Statistical Metabolomics Analysis – An R Tool. Analytical Chemistry 88, 6334 – 6341.                                                                                           
Yu-Jen Liang, Chih-Ting Yang, Chia-Wei Chen, Yin-Chun Lin, Shu-Yao Lin, Wen-Harn Pan, Yi Sheng Wang and Hsin-Chou Yang (2025/10). SMART 2.0 Statistical Metabolomics Analysis: An R Tool 2.0. Analytical Chemistry 97(46), 25453 – 25468. https://doi.org/10.1021/acs.analchem.5c03225

*Correspondence: Dr. Hsin-Chou Yang, Institute of Statistical Science, Academia Sinica. No. 128, Sec. 2, Academia Road, Nankang, Taipei 115, Taiwan. (Fax) 886-2-27886833; (Tel) 886-2-27875686; (E-mail) hsinchou@stat.sinica.edu.tw.
