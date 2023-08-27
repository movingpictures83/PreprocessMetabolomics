# Objective:
#   The purpose of the script is to transpose metabolomics data to feed into PLUMA pipeline
#       input:  columns are samples, rows are compound IDs
#       output: columns are compound IDs, rows are samples
#

import pandas as pd
import numpy as np
import os

from sklearn.impute import SimpleImputer


def get_mapping_dict(metadata_path):
    #metadata_path = "metabolites_names.txt"

    name_col = 1 #started from 0
    metabolite_col=-1

    map_dict = {}

    with open(metadata_path, 'r') as f:
        f.readline()
        for line in f.readlines():
            line = line.strip("\n")
            row = line.split("\t")
            id = row[metabolite_col]
            if id!="":
                map_dict[id] = id + "__" + row[name_col].replace(",","").strip('"').replace(" ", "").strip('"')
    return map_dict

def normalize(abundance_df, samples_col="", toIndex=True):
    if toIndex:
        abundance_df.index = abundance_df[samples_col]
        del abundance_df[samples_col]
    samples = list(abundance_df.index)
    for sample in samples:
        abundance_df = abundance_df.apply(lambda row: row/sum(row), axis=1)
    return abundance_df


import PyPluMA
import PyIO

class PreprocessMetabolomicsPlugin:
    def input(self, inputfile):
        self.parameters = PyIO.readParameters(inputfile)
        self.metabolomics_path = PyPluMA.prefix()+"/"+self.parameters["metabolomics_path"]
        self.metadata_path = PyPluMA.prefix()+"/"+self.parameters["metadata_path"]
        #metabolomics_path = "metabolon_HMDB.csv"
    def run(self):
        pass

    def output(self, outputfile):
       #out_path =  "metabolon_HMDB_filtered.csv"
       #out_norm_path = "metabolon_norm.csv"
       out_path = outputfile+"_HMDB_filtered.csv"
       out_norm_path = outputfile+"_norm.csv"
       df = pd.read_csv(self.metabolomics_path, dtype=str)
       df = df[~df["COMP ID"].isnull()]
       # Rename COMP ID
       name_dict = get_mapping_dict(self.metadata_path)
       df["COMP ID"] = df["COMP ID"].apply(lambda x: name_dict[x])
       df.index = df["COMP ID"]
       del df["COMP ID"]
       # Impute missing values with minimum
       columns = df.columns
       for col in columns:
          df[col] = df[col].astype(str)
          df[col] = df[col].apply(lambda x: x.replace(",",""))
       df = df.astype(float)
       for col in columns:
          df[col] = df[col].fillna(df[col].min())
       df = df.transpose()
       df.to_csv(out_path)
       # Normalize:
       norm_df = normalize(df, toIndex=False)
       norm_df.to_csv(out_norm_path)
