# DiGAS

A novel computational model for identifying genomic elements (from single exon to entire genomic regions)  that are linked to a given phenotypic condition (i.e. a disease) investigating for single nucleotide polymorphisms (SNP).

Genomic elements are investigated by means of a newly introduced descriptor of genomic information, the generalized allele spectrum. In contrast to allele spectrum, the novel descriptor takes into account the complete set of SNPs of a region at once. Thus, frequency is computed at region level rather than at SNP level. Then, we introduce the differential generalized allele spectrum, which captures the differences in frequency allele spectra between healthy and ill sets (control and case respectively).

Statistical significance of a gene is then evaluated by means of a differential analysis between the healthy and ill portions of the population.


## Contents

The srs folder contains all code used:

* **SNPsMapping.R**: generate the table containing the associatoin between each gene a its SNPs
* **Pvalues_ADNI.ipynb**: Generalized allele specturm and Pvalues calculation on ADNI cohort
* **Pvalues_PPMI_swedd.ipynb**: Generalized allele specturm and Pvalues calculation on PPMI cohort
* **Classification_ADNI.ipynb**: classification method on ADNI cohort
* **Classification_PPMI_SWEDD.ipynb**: classification method on PPMI cohort

## Data

To test this workflow we used just genotyping and demographic data avaiable at ADNI(http://adni.loni.usc.edu) and PPMI(http://ppmi-info.org) data portal.

## Usage

The project was developed using python inside jupyter notebook. For the mapping of SNPs to a gene was used R that can be lunched inside the notebook as well. 


## Contributing
Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

Please make sure to update tests as appropriate.

## Support

You can contact directly the authors by their e-mail addresses.

## Authors and acknowledgment

Guglielmo Cerri (cerriguglielmo@gmail.com), Antonino Aparo(antonino.aparo@univr.it) and Vincenzo Bonnici (vincenzo.bonnici@univr.it) created the workflow. 

## License
[MIT](https://choosealicense.com/licenses/mit/)

