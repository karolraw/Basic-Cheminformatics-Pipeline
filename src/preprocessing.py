from rdkit import Chem
from rdkit.Chem import SaltRemove, AllChem

def standardize_smiles(smiles: str) -> str:
  # Removing salts and neutralizing molecules in our dataset
  try:
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
      return None

    # remove salts\
    remover = SaltRemover.SaltRemover()
    mol = remover.StripMol(mol)

    # change all molecules to canonical SMILES
    return Chem.MolToSmiles(mol, canonical=True)

  except:
    return None

def preprocess_data(input_csv: str, output_csv: str):
  #Load the raw data, standardize SMILES, save clean data
  import pandas as pd

  df = pd.read_csv(input_csv)
  df["standardized_smiles"] = df["smiles"].apply(standardize_smiles)

  #Drop invalid rows
  df = df.dropna(subset=["standardized_smiles"])

  df.to_csv(output_csv, index=False)
  print(f"Data has been preprocessed and saved to {output.csv}")

if __name__ == "__main__":
  preprocess_data("data/raw/egfr_bioactivity_full.csv", "data/processed/egfr_cleaned.csv")
