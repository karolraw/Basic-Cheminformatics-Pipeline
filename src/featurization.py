from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors, AllChem
from rdkit.rdMolDescriptors import GetMorganFingerprintAsBitVect
import numpy as np
import pandas as pd
import os
import sys

project_root = os.path.dirname(os.path.abspath(__file__))
root_folder = os.path.abspath(os.path.join(project_root, '..'))
os.chdir(root_folder)

def compute_descriptors(smiles: str):
  # Compute descriptors from a given SMILES string
  mol = Chem.MolFromSmiles(smiles)
  if not mol:
    return None
  return {
    "MolWT": Descriptors.MolWt(mol),
    "TPSA": Descriptors.TPSA(mol),
    "NumHDonors": Descriptors.NumHDonors(mol),
    "NumHAcceptors": Descriptors.NumHAcceptors(mol),
    "LogP": Descriptors.MolLogP(mol)
  }

def compute_fingerprint(smiles: str, n_bits = 2048):
  # Compute Morgan fingerprints for a given SMILES
  mol = Chem.MolFromSmiles(smiles)
  if not mol:
    return None
  fp = GetMorganFingerprintAsBitVect(mol, radius=2, nBits=n_bits)
  return list(fp)

def featurize_data(input_csv: str, output_csv: str):
  # Input our data, compute descriptors and fingerprints, then save it
  df = pd.read_csv(input_csv)

  #Descriptors
  descriptors = df['canonical_smiles'].apply(compute_descriptors)
  descriptors_df = pd.DataFrame(descriptors.tolist())

  #Fingerprints
  fingerprints = df['canonical_smiles'].apply(compute_fingerprint)
  fingerprints_df = pd.DataFrame(fingerprints.tolist(), columns=[f"fp_{i}" for i in range(2048)])

  #Combine all the dataframes together
  final_df = pd.concat([df, descriptors_df, fingerprints_df], axis=1)

  final_df.to_csv(output_csv, index=False)
  print(f"Featurized data saved to {output_csv}")

if __name__ == "__main__":
  input_path = os.path.abspath("data/processed/egfr_cleaned.csv")
  output_path = os.path.abspath("data/processed/egfr_features.csv")  

  featurize_data(input_path, output_path)
