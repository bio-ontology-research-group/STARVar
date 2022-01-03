# STARVar

STARVar:Symptom based Tool for Automatic Ranking of Variants using evidence from literature and genomes

How to use STARVar?
To run STARVar, please follow the steps below:
0. Download the required input files from ....

1. An input VCF file first should be filtered for the candiate variants based on Allele Frequency and/or mode of inheritance.
We suggest users to use one of the efficient filtering tools slivar for the filtering step. Slivar is avaiulable from https://github.com/brentp/slivar

2. The filtered VCF file must be annotated by using the Variant Effect Predictor (VEP) tool with the required information for analyzing data efficiently. Follow the steps described https://github.com/bio-ontology-research-group/pavs_annotation to run the VEP as CWL.

3. run STARVar as follows:
format: python STARVar.py input.vcf option symptoms

option: STARVar combines scores from the literature based model and pathogenicity scores from either PolyPhen-2 or SIFT for variant ranking. Please choose one of them. Use s for SIFT and use p for PolyPhen-2 

symptoms: should be provided as a semi-colon separated single list of symptoms either in free text format or as HPO IDs. For example, "HP:XXXXXX;brain atrophy". Please dont forget to use quotation (") if you provide symtoms in free text format. 
