# DiGAS

A new methodology to prioritizes genes in relation to their capability of being biomarkers for a specific category, or for a group of two or more categories.

## Overview

Neurodegenerative diseases produce progressive loss of structure or function of neurons
and their causes can be found in the brain at many different levels of neuronal circuitry, ranging from molecular to systemic. Alzheimer's disease and Parkinson's disease are among the most common types of neurodegeneration. Evidences of genetic factors have been found for both diseases, however, because of their complexity their genetic diagnosis and treatment are still open challenges.
Alzheimer's disease is the most common cause of dementia. Nowadays, the classification for this disease is very challenging and have been presented 
Different algorithms have been proposed for investigating the genetic causes of neurodegenerative diseases. Most of them use single data modality or deep learning (DL) to solve the problem. Here, using Alzheimer’s disease neuroimaging initiative (ADNI) dataset, 
we presented a Differential generalized gene allele spectrum as descriptor of the genetic data together with a computational model that exploits it for recognizing affected individuals in genetic studies. 
The model is able to predict, with high accuracy, the onset of Alzheimer's disease and Parkinson's disease using only single-nucleotide polymorphisms.

Here is showed a method illustration


![Image](doc/method_illustration.PNG)


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

