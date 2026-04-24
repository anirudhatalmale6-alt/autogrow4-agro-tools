#!/bin/bash
# Setup script for agricultural_reactions library in AutoGrow4
# Run from ~/final_project/

set -e

AUTOGROW_DIR="${HOME}/final_project/autogrow4"
RXN_LIB_DIR="${AUTOGROW_DIR}/autogrow/operators/mutation/smiles_click_chem/reaction_libraries"
AGRO_DIR="${RXN_LIB_DIR}/agro_rxns"
ALL_RXNS_DIR="${RXN_LIB_DIR}/all_rxns"

echo "Setting up agricultural reactions library..."
echo "  AutoGrow4: ${AUTOGROW_DIR}"
echo "  Reaction libraries: ${RXN_LIB_DIR}"

# Create agro_rxns directory
mkdir -p "${AGRO_DIR}/complementary_mol_dir"

# Copy JSON files
cp Agro_Rxns_rxn_library.json "${AGRO_DIR}/"
cp Agro_Rxns_functional_groups.json "${AGRO_DIR}/"

# Copy new agricultural complementary molecule files
for f in complementary_mols/*.smi; do
    cp "$f" "${AGRO_DIR}/complementary_mol_dir/"
done

# Copy existing click-chem and robust .smi files from all_rxns
if [ -d "${ALL_RXNS_DIR}/complementary_mol_dir" ]; then
    for f in "${ALL_RXNS_DIR}/complementary_mol_dir"/*.smi; do
        basename=$(basename "$f")
        if [ ! -f "${AGRO_DIR}/complementary_mol_dir/${basename}" ]; then
            cp "$f" "${AGRO_DIR}/complementary_mol_dir/"
        fi
    done
    echo "  Copied existing .smi files from all_rxns"
else
    echo "  WARNING: all_rxns complementary_mol_dir not found at:"
    echo "    ${ALL_RXNS_DIR}/complementary_mol_dir"
    echo "  You may need to copy .smi files manually"
fi

# Count files
n_json=$(ls "${AGRO_DIR}"/*.json 2>/dev/null | wc -l)
n_smi=$(ls "${AGRO_DIR}/complementary_mol_dir"/*.smi 2>/dev/null | wc -l)

echo ""
echo "Done! Agricultural reactions library installed:"
echo "  ${AGRO_DIR}/"
echo "  JSON files: ${n_json}"
echo "  .smi files: ${n_smi}"
echo ""
echo "To use in AutoGrow4, set:"
echo '  --rxn_library custom'
echo "  --rxn_library_file ${AGRO_DIR}/Agro_Rxns_rxn_library.json"
echo "  --function_group_library ${AGRO_DIR}/Agro_Rxns_functional_groups.json"
echo "  --complementary_mol_directory ${AGRO_DIR}/complementary_mol_dir/"
