from rcsbapi.search import AttributeQuery, search_attributes as attrs
import json

# Load retrieved EGFR PDB IDs
with open("egfr_pdb_ligands.json", "r") as f:
    pdb_ids = json.load(f)

# List of approved EGFR inhibitors
approved_drugs = {
    "Gefitinib": "DB00530",
    "Erlotinib": "DB00530",
    "Osimertinib": "DB09330",
    "Lapatinib": "DB01229",
    "Afatinib": "DB08916",
    "Dacomitinib": "DB11763"
}

# Query to get ligand details
filtered_pdbs = []

for pdb_id in pdb_ids:
    query = AttributeQuery(
        attribute="rcsb_entry_container_identifiers.entry_id",
        operator="exact_match",
        value=pdb_id
    )

    # Fetch results
    result = list(query())

    for entry in result:
        ligands = entry.get("nonpolymer_entities", [])
        for ligand in ligands:
            ligand_name = ligand.get("chem_comp.name", "").lower()
            if any(drug.lower() in ligand_name for drug in approved_drugs.keys()):
                filtered_pdbs.append(pdb_id)
                break  # Stop checking once a match is found

# Save the filtered list
with open("egfr_pdb_filtered.json", "w") as f:
    json.dump(filtered_pdbs, f, indent=4)

print(f"Filtered {len(filtered_pdbs)} PDB structures with marketed EGFR inhibitors.")
