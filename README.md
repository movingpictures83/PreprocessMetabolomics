# PreprocessMetabolomics
# Language: Python
# Input: TXT
# Output: PREFIX
# Tested with: PluMA 1.1, Python 3.6
# Dependency: pandas==1.1.5

Metabolomics data preprocessor.

Input is a TXT file of keyword-value pairs (Tab-delimited):

metabolomics_path: Metabolite concentrations, CSV
metadata_path: Metadata,TXT

Output metabolomics data is in unnormalized and normalized CSV format, using the user-specified prefix.
