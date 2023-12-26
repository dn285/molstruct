from flask import Flask, request, jsonify
from flask_cors import CORS

#from openbabel import pybel
#import pubchempy as pcp

from collections import Counter

app = Flask(__name__)
CORS(app, resources={r"/api/*": {"origins": "https://molstruct-c555.vercel.app"}})

@app.route('/')
def home():
    return "Molstruct backend up and running!"

def get_molecular(atoms):
    d_atoms = Counter(list(zip(*atoms))[1])

    formula = ""
    if 'C' in d_atoms:
        formula += f"C{d_atoms['C']}" if d_atoms['C'] > 1 else 'C'
        del d_atoms['C']
    if 'H' in d_atoms:
        formula += f"H{d_atoms['H']}" if d_atoms['H'] > 1 else 'H'
        del d_atoms['H']
    for atm in sorted(d_atoms.keys()):
        formula += f"{atm}{d_atoms[atm]}" if d_atoms[atm] > 1 else atm
    return formula

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
        return jsonify({
            "molecular": "TEST",
            "smiles": "TEST",
            "stdinchi": "TEST",
            "stdinchikey": "TEST",
            "iupac": "TEST",
        })

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
        my_molecule = pybel.readstring('sdf', molfile)
        smiles = my_molecule.write('smi')

        # Obtain IUPAC name from PubChem
        queried_mol = pcp.get_compounds(smiles, 'smiles')[0]

        return jsonify({
            "molecular": queried_mol.molecular_formula,
            "smiles": smiles,
            "stdinchi": my_molecule.write('inchi'),
            "stdinchikey": my_molecule.write('inchikey'),
            "iupac": queried_mol.iupac_name
        })
    
    except Exception as e:
        print(f"Backend error: {e}")
        return jsonify({"error": "Error processing data"}), 500

if __name__ == '__main__':
    app.run(debug=True)