from flask import Flask, request, jsonify
from flask_cors import CORS

app = Flask(__name__)
CORS(app)

@app.route('/')
def home():
    return "Molstruct backend up and running!"

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

        print(data['atoms'])

        atoms = [(atom['id'], atom['type']) for atom in data['atoms']]
        bonds = [(bond['start']['id'], bond['end']['id']) for bond in data['bonds']]

        # TODO: Code to compute SMILES and others
        smiles = get_smiles(atoms, bonds)
        smarts = get_smarts(atoms, bonds)
        stdinchi = get_stdinchi(atoms, bonds)
        stdinchikey = get_stdinchikey(atoms, bonds)
        iupac = get_iupac(atoms, bonds)

        return jsonify({
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