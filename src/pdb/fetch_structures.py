from rcsbapi.search import AttributeQuery
from rcsbapi.search import search_attributes as attrs

# Search for EGFR structures by UniProt Accession ID
query_egfr = AttributeQuery(
    attribute="rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession",
    operator="exact_match",
    value="P00533"
)

# Filter for structures that contain bound ligands (chemical components)
query_ligands = AttributeQuery(
    attribute="rcsb_entry_container_identifiers.entry_id",
    operator="exists"
)

# Combine queries: EGFR structures that also have bound ligands
query_combined = query_egfr & query_ligands

# Fetch results
pdb_ids = list(query_combined())

print(f"Found {len(pdb_ids)} EGFR PDB structures with bound ligands.")

# Save results
import json

with open("egfr_pdb_ligands.json", "w") as f:
    json.dump(pdb_ids, f, indent=4)

print("Saved filtered PDB IDs to egfr_pdb_ligands.json")
