from flask import Flask, request, jsonify
from flask_cors import CORS

from rdkit import Chem
import pubchempy as pcp

from collections import Counter

app = Flask(__name__)
CORS(app, resources={r"/convert": {"origins": "https://molstruct-c555.vercel.app"}})

@app.route('/')
def home():
    return "Molstruct backend up and running!"

def get_molfile(atoms, bonds):
    lines = [
        '', 
        'Molstruct by Dan Ni', 
        'Visit https://github.com/dn285/molstruct for more information'
    ]

    # Counts line
    counts_line = [len(atoms), len(bonds)] + [0] * 7
    counts_line = ''.join(f"{field:>3}" for field in counts_line) 
    counts_line += "  0999" + " V2000"
    lines.append(counts_line)
    
    # Atoms block
    for atom in atoms:
        atom_line = [atom['x'], atom['y'], 0.0]
        atom_line = ''.join(f"{field:>10.4f}" for field in atom_line)
        atom_line += f" {atom['type']:>3}" + "  0" # Edit to implement charge
        atom_line += "  0" * 11
        lines.append(atom_line)

    # Bonds block
    for start, end, multiplicity in bonds:
        bond_line = [
            start+1, 
            end+1, 
            multiplicity, # Edit to implement multiple bonds
            0 # Edit to implement stereochemistry
        ]
        bond_line += [0] * 3
        bond_line = ''.join(f"{field:3}" for field in bond_line)
        lines.append(bond_line)

    # Terminator 
    lines.append("M  END")
    
    return '\n'.join(lines)

@app.route('/convert', methods=['POST'])
def compute_structural_data():
    try:
        data = request.json

        atoms = [(atom['id'], atom['type']) for atom in data['atoms']]

        # Remove duplicate bonds
        bonds = []
        for bond in data['bonds']:
            start, end = bond['start']['id'], bond['end']['id']
            multiplicity = bond['multiplicity']
            if (start, end) in bonds:
                continue
            if (end, start) in bonds:
                continue
            bonds.append((start, end, multiplicity))

        # Compute SMILES and others
        molfile = get_molfile(data['atoms'], bonds)
        my_molecule = Chem.MolFromMolBlock(molfile) 
        smiles = Chem.MolToSmiles(my_molecule, isomericSmiles=False)
        smarts = Chem.MolToSmarts(my_molecule, isomericSmiles=False)

        print(smiles)

        # Obtain IUPAC name from PubChem
        queried_mol = pcp.get_compounds(smiles, 'smiles')[0]

        return jsonify({
            "molecular": queried_mol.molecular_formula,
            "smiles": smiles,
            "smarts": smarts,
            "iupac": queried_mol.iupac_name
        })
    
    except Exception as e:
        print(f"Backend error: {e}")
        return jsonify({"error": "Error processing data"}), 500

if __name__ == '__main__':
    app.run(debug=True)