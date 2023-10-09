# Description
Pipelines and scripts used to perform the bioinformatics analyses examining the parental factors driving pathogenicity of variably expressive variants.

# Repository Structure
Preprocessing of files are very similar for all cohorts, but scripts for each cohort are provided in 0_preprocessing.
Analyses of each cohort are primarily organized by the cohort of interest and files related to analyses for each cohort can be found in their respective directories.
A brief description of each repository is provided below:

**0_preprocessing**: This directory holds all of the data pre-processing scripts and is broken into separate directories for the three main types of pre-processed data used in the manuscript.  
  *1_variant_calling*: Scripts related to SNV calling in each cohort.  
  *2_kinship_PCs*: Scripts for calculating kinship coefficients between spouses and genetic principal components for each sample in each cohort.  
  *3_ROH*: Scripts for calculting runs of homozygosity in probands in SPARK and SSC.  
**1_SPARK**: This directory holds scripts related to analyses involving samples from the SPARK cohort, including those shown in Fig. 1A, 1C, 2C, and 3B  
**2_SSC**: This directory holds scripts related to analyses involving samples from the SSC cohort, including those shown in Fig. 1B, 2E-F, and 5A  
**3_16p12.1**: This directory holds scripts related to analyses involving samples from the Girirajan Lab 16p12.1 deletion cohort, including those shown in Fig. 2A-B and S3  
**4_UK Biobank**: This directory holds scripts related to analyses involving samples from the UK Biobank cohort, including those shown in Fig. 1D  
**5_MultiCohort**:This directory holds scripts related to analyses involving samples from multiple cohorts, including those shown in Fig. 4, S4, and S5  
**6_GeneDx**: This directory holds scripts related to analyses involving samples from the GeneDx cohort, including those shown in Fig. 5B-C  
**7_Modeling**:This directory holds scripts related to the simulation analysis shown in Fig. 3C-D and S2  

# Copyright/License
The code in this repository is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This code is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with the code.  If not, see <https://www.gnu.org/licenses/>.

# Contact
For questions or comments, please contact Corrine Smolen (ces6136@psu.edu) or Santhosh Girirajan (sxg47@psu.edu).
