# Molstruct

## About

Molstruct is a web application where users draw molecular structures and obtain the corresponding formula, structural data and IUPAC name.

## Getting Started

In the project directory: Run `npm start` to open the React frontend. If nothing opens, you may need to manually navigate to http://localhost:3000 (or the port you have set for React) in the web browser of your choice. If you are only using Molstruct for the basic canvas functionality, you may stop here.

To spin up the backend, which is required for accessing naming functionality: Navigate to the backend folder, install dependencies with `pip install -r requirements.txt`, then run the backend file with `python app.py`. It is recommended that this be done within a virtual environment.

## How it Works

The Molstruct backend constructs a MOL file that is translated into SMILES, StdInChI and StdInChIKey formats using OpenBabel. The molecular formula and IUPAC name are obtained by querying PubChem using the PubChemPy API with the aforementioned SMILES formula. As a result, IUPAC names are not guaranteed to be correct. In particular, if a molecule is not in the PubChem database, the displayed IUPAC name will be incorrect.

## Learn More

This project was bootstrapped with [Create React App](https://github.com/facebook/create-react-app). You can learn more in the [Create React App documentation](https://facebook.github.io/create-react-app/docs/getting-started).

To learn React, check out the [React documentation](https://reactjs.org/).

## License

This project is licensed under the GPL-3.0 license - see the [LICENSE](LICENSE) file for details.
