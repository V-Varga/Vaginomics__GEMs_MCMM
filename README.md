# Vaginomics__GEMs_MCMM

Author: Vi Varga

Project Period: August 2025

Course Code: FKBT190


## Repository Description

Workflow for creation of Genome-Scale Metabolic Models and a Microbial Community Metabolic Model to model the human vaginal environment


## Background & Project Description

Genome Scale Metabolic Models (GEMs) are computational reconstructions of an organism’s metabolism. While the curation of any GEM remains a massive undertaking, involving much literature review and experimental validation; in the last decade it has become increasingly straightforward to reconstruct a draft GEM for an organism based on DNA sequencing data. The improving quality and accessibility of GEMs has led to efforts to construct community-scale metabolic models (referred to as Microbial Community Metabolic Models [MCMMs]) by modeling metabolic exchange interactions between different GEMs. For this project I use MATLAB and Jupyter to create a draft human vaginal MCMM. 


## Workflow

The `conda` environment used to run all Jupyter Notebook-based code is included in this repository as the `env-GEM-MCMM.yml` file. 

It can be used to create a `conda` environment like so: 

```bash
conda env create -f env-GEM-MCMM.yml
```

Protein FASTA files downloaded from the NCBI were cleaned in preparation for use in GEM creation with the `assignFASTAheaders_GEMs.py` script.

GEMs were created from using the [RAVEN Toolbox](https://github.com/SysBioChalmers/RAVEN) in MATLAB from RefSeq protein FASTAs downloaded from the NCBI:
 - *Lactobacillus crispatus*
 - *Lactobacillus iners*
 - *Prevotella bivia*
 - *Gardnerella vaginalis*
 - *Trichomonas vaginalis*

A KEGG version of [Human-GEM](https://github.com/SysBioChalmers/Human-GEM) was also created in MATLAB. The original Human-GEM XML file I downloaded from the GitHub repository is included in the `Data/` directory for transparency.

The MCMM was created using [PyCoMo](https://github.com/univieCUBE/PyCoMo) in a Jupyter Notebook. 

All created and modified SBML-formatted XML files of the GEMs and the MCMM are included in the `Results_MCMMs/` directory. All scripts needed to run the MATLAB Live script and Jupyter Notebook are included in the `Scripts/` directory. All data files necessary to run the MATLAB Live script are included in the `Data/` directory.


## References

General methodology & background: 
 - Ardalani, O., Phaneuf, P. v., Mohite, O. S., Nielsen, L. K., & Palsson, B. O. (2024).  Pangenome reconstruction of Lactobacillaceae metabolism predicts species-specific metabolic traits. MSystems, 9(7). https://doi.org/10.1128/MSYSTEMS.00156-24/SUPPL_FILE/MSYSTEMS.00156-24-S0008.DOCX
 - Diener, C., & Gibbons, S. M. (2023). More is Different: Metabolic Modeling of Diverse Microbial Communities. MSystems, 8(2). https://doi.org/10.1128/msystems.01270-22
 - Dillard, L. R., Glass, E. M., Lewis, A. L., Thomas-White, K., Papin, J. A., & Bucci, V. (2025). Metabolic Network Models of the Gardnerella Pangenome Identify Key Interactions with the Vaginal Environment. https://doi.org/10.1128/msystems.00689-22
 - Fondi, M., & Liò, P. (2015). Genome-scale metabolic network reconstruction. Methods in Molecular Biology, 1231, 233–256. https://doi.org/10.1007/978-1-4939-1720-4_15
 - Lambert, A., Budinich, M., Mahé, M., Chaffron, S., & Eveillard, D. (2024). Community metabolic modeling of host-microbiota interactions through multi-objective optimization. IScience, 27(6). https://doi.org/10.1016/j.isci.2024.110092
 - Mendoza, S. N., Olivier, B. G., Molenaar, D., & Teusink, B. (2019). A systematic assessment of current genome-scale metabolic reconstruction tools. Genome Biology, 20(1). https://doi.org/10.1186/s13059-019-1769-1
 - Mengoni, A., Galardini, M., & Fondi, M. (n.d.). Bacterial Pangenomics Methods and Protocols Methods in Molecular Biology 1231. http://www.springer.com/series/7651
 - Predl, M., Mießkes, M., Rattei, T., & Zanghellini, J. (2024). PyCoMo: a python package for community metabolic model creation and analysis. Bioinformatics, 40(4). https://doi.org/10.1093/bioinformatics/btae153
 - Quinn-Bohmann, N., Carr, A. v., Diener, C., & Gibbons, S. M. (2025). Moving from genome-scale to community-scale metabolic models for the human gut microbiome. Nature Microbiology, 10(5), 1055–1066. https://doi.org/10.1038/s41564-025-01972-2
 - Robinson, J. L., Kocabaş, P., Wang, H., Cholley, P. E., Cook, D., Nilsson, A., Anton, M., Ferreira, R., Domenzain, I., Billa, V., Limeta, A., Hedin, A., Gustafsson, J., Kerkhoven, E. J., Svensson, L. T., Palsson, B. O., Mardinoglu, A., Hansson, L., Uhlén, M., & Nielsen, J. (2020). An Atlas of Human Metabolism. Science Signaling, 13(624), eaaz1482. https://doi.org/10.1126/SCISIGNAL.AAZ1482
 - Wang, H. I., Marciš auskas, S. I., Sá nchez ID, B. J., Domenzain, I. I., Hermansson, D., Agren, R., NielsenID, J., & KerkhovenID, E. J. (2018). RAVEN 2.0: A versatile toolbox for metabolic network reconstruction and a case study on Streptomyces coelicolor. https://doi.org/10.1371/journal.pcbi.1006541

Data sources: 
 - Human-GEM: https://github.com/SysBioChalmers/Human-GEM
 - Organism FASTA files were downloaded from the NCBI

Note that ChatGPT was used to aid in coding, primarily in the writing of MATLAB functions. 
