# DiGAS

A novel computational model for identifying genomic elements (from single exon to entire genomic regions)  that are linked to a given phenotypic condition (i.e. a disease) investigating for single nucleotide polymorphisms (SNP).

Genomic elements are investigated by means of a newly introduced descriptor of genomic information, the generalized allele spectrum. In contrast to allele spectrum, the novel descriptor takes into account the complete set of SNPs of a region at once. Thus, frequency is computed at region level rather than at SNP level. Then, we introduce the Differential Generalized Allele Spectrum (DiGAS), which captures the differences in frequency allele spectra between healthy and ill sets (control and case respectively). Statistical significance of a gene is then evaluated by means of a differential analysis between the healthy and ill portions of the population.

DiGAS is distributed as a command-line tool developed in Python and it is available for both Windows and Unix systems. It uses PLINK to retrieve data about genotyping and takes as input the coordinates of genomic regions to be analyzed.


## Installation and usage

DiGAS can run on Windows or Unix system where Python 3 (or higher) was previously installed.  
DiGAS depends on PLINK was previously installed , so we first need to install this package through the remotes package (in this way we can install a fixed version of Seurat). After Seurat installation, we can install Stardust from our GitHub repository. 
Thus, the user will need Cython to be installed to correctly build GRAFIMO. Cython can be obtained via pip


DiGAS depends on a number of Python packages. To install all the required packages:


## Data

To test DiGAS we used just genotyping and demographic data avaiable at ADNI(http://adni.loni.usc.edu) data portal.

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

