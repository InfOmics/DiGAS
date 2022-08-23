# DiGAS

A novel computational model for identifying genomic elements (from single exon to entire genomic regions)  that are linked to a given phenotypic condition (i.e. a disease) investigating for single nucleotide polymorphisms (SNP).

Genomic elements are investigated by means of a newly introduced descriptor of genomic information, the generalized allele spectrum. In contrast to allele spectrum, the novel descriptor takes into account the complete set of SNPs of a region at once. Thus, frequency is computed at region level rather than at SNP level. Then, we introduce the Differential Generalized Allele Spectrum (DiGAS), which captures the differences in frequency allele spectra between healthy and ill sets (control and case respectively). Statistical significance of a gene is then evaluated by means of a differential analysis between the healthy and ill portions of the population.

DiGAS is distributed as a command-line tool developed in Python and it is available for both Windows and Unix systems. It uses PLINK to retrieve data about genotyping and takes as input the coordinates of genomic regions to be analyzed.


## Installation and usage

DiGAS can be easly downloaded from GitHub `src` directory (Windows or Unix version). PLINK tool is also contained in both directory.

### System requirements

DiGAS can run on Windows or Unix system where Python 3 (or higher) was previously installed.

DiGAS depends on PLINK. PLINK executable must be contained in the same folder as the DiGAS executable.
A PLINK version is contained in the `src` directory for both Windows and Unix versions. Alternatively you can download it from the official website: https://www.cog-genomics.org/plink/.

DiGAS depends on a number of Python packages. To install all the required packages:

```
pip3 install pandas
pip3 install numpy
pip3 install argparse
pip3 install csv
pip3 install subprocess
pip3 install shutil
pip3 install re
```

### Usage
TO DO 
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

