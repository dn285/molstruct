from flask import Flask, request, jsonify
from flask_cors import CORS

from collections import Counter

app = Flask(__name__)
CORS(app)

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

def get_smiles(atoms, bonds):
    print(atoms)
    print(bonds)

    return "Coming soon!"

def get_smarts(atoms, bonds):
    return "Coming soon!"

def get_stdinchi(atoms, bonds):
    return "Coming soon!"

def get_stdinchikey(atoms, bonds):
    return "Coming soon!"

def get_iupac(atoms, bonds):
    return "Coming soon!"

@app.route('/convert', methods=['POST'])
def compute_structural_data():
    try:
        data = request.json

        atoms = [(atom['id'], atom['type']) for atom in data['atoms']]
        bonds = [(bond['start']['id'], bond['end']['id']) for bond in data['bonds']]

        # TODO: Code to compute SMILES and others
        molecular = get_molecular(atoms)
        smiles = get_smiles(atoms, bonds)
        smarts = get_smarts(atoms, bonds)
        stdinchi = get_stdinchi(atoms, bonds)
        stdinchikey = get_stdinchikey(atoms, bonds)
        iupac = get_iupac(atoms, bonds)

        return jsonify({
            "molecular": molecular,
            "smiles": smiles,
            "smarts": smarts,
            "stdinchi": stdinchi,
            "stdinchikey": stdinchikey,
            "iupac": iupac
        })
    except Exception as e:
        print(f"Backend error: {e}")
        return jsonify({"error": "Error processing data"}), 500

if __name__ == '__main__':
    app.run(debug=True)