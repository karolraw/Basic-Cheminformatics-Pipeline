from rdkit import Chem
from rdkit.Chem import AllChem
import os
import sys

project_root = os.path.dirname(os.path.abspath(__file__))
os.chdir(project_root)

def preprocess_data(input_csv: str, output_csv: str):
  #Load the bioactivity data, do some minor processing, and add an activity column for active or inactive
  import pandas as pd
  print(f"Reading from: {input_csv}")

  if not os.path.exists(input_csv):
    raise FileNotFoundError(f"Input file not found: {input_csv}")  
  
  df = pd.read_csv(input_csv)

  #Drop rows with missing values
  df = df.dropna()

  df['standard_value'] = pd.to_numeric(df['standard_value'], errors='coerce')
  threshold=1000
  df['bioactivity_class'] = df['standard_value'].apply(lambda x: 'active' if x < threshold else 'inactive')

  df.to_csv(output_csv, index=False)
  print(f"Data has been preprocessed and saved to {output_csv}")

if __name__ == "__main__":
  input_csv = os.path.abspath("../data/raw/egfr_bioactivity_full.csv")
  output_csv = os.path.abspath("../data/processed/egfr_cleaned.csv")
  preprocess_data(input_csv, output_csv)
