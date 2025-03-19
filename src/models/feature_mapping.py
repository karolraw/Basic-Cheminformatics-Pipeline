import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw, AllChem

#Read our previously generated data
data = pd.read_csv('data/processed/egfr_features.csv')
feature_importances = pd.read_csv('data/processed/feature_importances.csv')

#Take the top 20 features
top_features = feature_importances['Feature'].tolist()[:20]
top_bit_indices = [int(feature.split(' ')[-1]) for feature in top_features]

#Load the first molecule in our saved data
smiles_to_analyze = data['canonical_smiles'][0]
mol_to_analyze = Chem.MolFromSmiles(smiles_to_analyze)

#Mapping the bit info to substructure
bit_info = {}
fp = AllChem.GetMorganFingerprintAsBitVect(mol_to_analyze, radius=2, nBits=2048, bitInfo=bit_info)

substructure_mols = []
bit_legends = []

for bit in top_bit_indices:
  if bit in bit_info:
    atom_ids, radius = bit_info[bit][0]
    env = Chem.FindAtomEnvironmentOfRadiusN(mol_to_analyze, radius, atom_ids)
    amap = {}
    submol = Chem.PathToSubmol(mol_to_analyze, env, atomMap=amap)
    substructure_mols.append(submol)
    bit_legends.append(f"Bit {bit}")

#Visualize and save the substructures
img = Draw.MolsToGridImage(substructure_mols, molsPerRow=5, subImgSize=(300,300), legends=bit_legends)
img.save('results/top_20_substructures.png')
