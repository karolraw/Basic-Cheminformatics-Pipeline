import os
import sys
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem, Draw, DataStructs
import seaborn as sns
import matplotlib.pyplot as plt

#Make sure it runs from the root directory
project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "../../"))
os.chdir(project_root)

#Load data
data = pd.read_csv("data/processed/egfr_features.csv")

#Select a sample of active/inactive molecules
active_mols = data[data["bioactivity_class"] == "active"]["canonical_smiles"].head(5)
inactive_mols = data[data["bioactivity_class"] == "inactive"]["canonical_smiles"].head(5)

#Convert SMILES to RDKit objects
active_mols = [Chem.MolFromSmiles(smiles) for smiles in active_mols]
inactive_mols = [Chem.MolFromSmiles(smiles) for smiles in inactive_mols]

#Visualize the compounds
Draw.MolsToGridImage(active_mols, molsPerRow=5, subImgSize=(300,300), legends=["Active"]*5).save("src/models/active_molecules.png")
Draw.MolsToGridImage(inactive_mols, molsPerRow=5, subImgSize=(300,300), legends=["Inactive"]*5).save("src/models/inactive_molecules.png") 

#Load feature importances
feature_importances = pd.read_csv("data/processed/feature_importances.csv")
bit_indices = feature_importances["Feature"].str.extract(r"(\d+)")[0].astype(int).tolist()

smiles_list = data["canonical_smiles"].tolist()
mols = [Chem.MolFromSmiles(smiles) for smiles in smiles_list if smiles is not None]

bit_info = {}
for mol in mols[:1]:
  if mol is not None:
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=2048, bitInfo = bit_info)

bit_info_filtered = {bit: bit_info[bit] for bit in bit_indices if bit in bit_info}

substructure_dir = "data/processed/substructures/"
os.makedirs(substructure_dir, exist_ok=True)

for bit, info in bit_info_filtered.items():
  atom_id, radius = info[0]
  env = Chem.FindAtomEnvironmentOfRadiusN(mols[0], radius, atom_id)
  submol = Chem.PathToSubmol(mols[0], env)
  img = Draw.MolToImage(submol, size=(200,200))
  img.save(f"{substructure_dir}substructure_{bit}.png")

active_mols = [Chem.MolFromSmiles(smiles) for smiles in data[data["bioactivity_class"]=="active"]["canonical_smiles"]]
inactive_mols = [Chem.MolFromSmiles(smiles) for smiles in data[data["bioactivity_class"] == "inactive"]["canonical_smiles"]]

#Count the substructure occurences
def count_substructures(mols, bit_info_filtered):
  substructure_counts = {bit: 0 for bit in bit_indices}
  for mol in mols:
    if mol is not None:
      fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=2048)
      for bit in bit_indices:
        if fp[bit]:
          substructure_counts[bit] += 1
  return substructure_counts

active_counts = count_substructures(active_mols, bit_info_filtered)
inactive_counts = count_substructures(inactive_mols, bit_info_filtered)

#Create a new dataframe for comparison
substructure_df = pd.DataFrame({
  "Bit Index": bit_indices,
  "Active Count": [active_counts[bit] for bit in bit_indices],
  "Inactive Count": [inactive_counts[bit] for bit in bit_indices]
})
substructure_df.to_csv("data/processed/substructure_comparison_all.csv", index=False)

#Look through the top 50
top_n = 50
substructure_df_top = substructure_df.sort_values(by="Active Count", ascending=False).head(top_n)

#create a heatmap
plt.figure(figsize=(12,8))
sns.heatmap(substructure_df_top[["Active Count", "Inactive Count"]].T, annot=True, fmt="d", cmap="coolwarm", xticklabels=substructure_df_top["Bit Index"])
plt.title(f"Substructure Frequency in Active vs Inactive Compounds (Top {top_n})")
plt.xlabel("Bit Index")
plt.ylabel("Count")
plt.savefig("data/processed/substructure_heatmap_all.png")
plt.show()
