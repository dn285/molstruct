import React, { useState, useRef, useEffect } from 'react';
import './App.css';

function App() {
    const elements = ['H', 'B', 'C', 'N', 'O', 'P', 'S', 'F', 'Cl', 'Br', 'I'];

    const [selectedTool, setSelectedTool] = useState('C');
    const [atoms, setAtoms] = useState([]);
    const [bonds, setBonds] = useState([]);
    const [selectedAtoms, setSelectedAtoms] = useState([]);
    const canvasRef = useRef(null);
    const atomRadius = 8;

    //const [adjacencyMatrix, setAdjacencyMatrix] = useState([]);
    //const [atomData, setAtomData] = useState([]);
    const [structuralData, setStructuralData] = useState({
        molecular: '',
        smiles: '',
        stdinchi: '',
        stdinchikey: '',
        iupac: 'Coming soon!'
    });

    const clearCanvas = () => {
        setAtoms([]);
        setBonds([]);
    };

    useEffect(() => {
        const canvas = canvasRef.current;
        const context = canvas.getContext('2d');

        // Clear canvas
        context.clearRect(0, 0, canvas.width, canvas.height);

        // Draw atoms
        atoms.forEach(atom => {
            const x = atom.x
            const y = atom.y

            context.font = "16px Arial";
            context.textAlign = "center";
            context.textBaseline = "middle";
            context.fillText(atom.type, x, y);
        });

        // Draw bonds
        bonds.forEach(bond => {
            context.beginPath();
            context.moveTo(bond.start.x, bond.start.y);
            context.lineTo(bond.end.x, bond.end.y);
            context.stroke();
        });

        updateStructuralData();
    }, [atoms, bonds]);

    // Handlers

    const handleCanvasClick = (event) => {
        const canvas = canvasRef.current;
        const rect = canvas.getBoundingClientRect();
        const scaleX = canvas.width / rect.width;
        const scaleY = canvas.height / rect.height;

        const x = (event.clientX - rect.left) * scaleX;
        const y = (event.clientY - rect.top) * scaleY;

        const clickedAtom = atoms.find(atom =>
            Math.abs(atom.x - x) < 10 && Math.abs(atom.y - y) < 10);

        if (elements.includes(selectedTool)) {
            addAtom(selectedTool, x, y);
            setSelectedAtoms([]);
        }
        else if (selectedTool === 'bond') {
            if (clickedAtom) {
                if (selectedAtoms.length < 1) {
                    setSelectedAtoms([clickedAtom]);
                }
                else {
                    if (selectedAtoms[0].id !== clickedAtom.id) {
                        const newBond = createBond(selectedAtoms[0], clickedAtom, atomRadius);
                        setBonds([...bonds, newBond])
                    }

                    setSelectedAtoms([]);
                }
            }
        }
    };

    const addAtom = (atomType, x, y) => {
        const newAtom = {
            id: atoms.length,
            type: atomType,
            x: x,
            y: y
        };
        setAtoms([...atoms, newAtom]);
    };

    const createBond = (startAtom, endAtom, atomRadius) => {
        const dx = endAtom.x - startAtom.x;
        const dy = endAtom.y - startAtom.y;
        const dist = Math.sqrt(dx * dx + dy * dy);

        const vx = dx / dist;
        const vy = dy / dist;
        const offsetX = vx * atomRadius;
        const offsetY = vy * atomRadius;

        return {
            start: {
                id: startAtom.id,
                x: startAtom.x + offsetX,
                y: startAtom.y + offsetY
            },
            end: {
                id: endAtom.id,
                x: endAtom.x - offsetX,
                y: endAtom.y - offsetY
            }
        };
    };

    const updateStructuralData = async () => {
        try {
            const data = await fetchStructuralData();
            if (data) {
                setStructuralData(data);
            }
        }
        catch (error) {
            console.error("Error during conversion:", error);
        }
    };

    const fetchStructuralData = async () => {
        try {
            const response = await fetch('http://127.0.0.1:5000/convert', {
                method: 'POST',
                headers: {
                    'Content-Type': 'application/json',
                },
                body: JSON.stringify({ atoms, bonds })
            });
            const structuralData = await response.json();
            return structuralData;
        }
        catch (error) {
            console.error("Error during fetching:", error);
        }
    };

    // The actual page

    return (
        <div className="App">
            <div className="top-container">
                <div className="sidebar">
                    <h1>Welcome to molstruct!</h1>
                    <p>Draw your molecule onto the canvas on the right and see the corresponding molecular formula below.</p>
                    <div className="output-names">
                        <h2>Outputs</h2>
                        <div className="name-group">
                            <label>Molecular Formula</label>
                            <input type="text" readOnly value={structuralData.molecular} />
                        </div>
                        <div className="name-group">
                            <label>SMILES</label>
                            <input type="text" readOnly value={structuralData.smiles} />
                        </div>
                        <div className="name-group">
                            <label>StdInChI</label>
                            <input type="text" readOnly value={structuralData.stdinchi} />
                        </div>
                        <div className="name-group">
                            <label>StdInChIKey</label>
                            <input type="text" readOnly value={structuralData.stdinchikey} />
                        </div>
                        <div className="name-group">
                            <label>IUPAC Name</label>
                            <input type="text" readOnly value={structuralData.iupac} />
                        </div>
                    </div>
                    <div>
                        <h2>Export</h2>
                        <p>Coming soon!</p>
                    </div>
                </div>
                <div className="canvas-container">
                    <div className="toolbar">
                        {elements.map(element => (
                            <button
                                key={element}
                                className={selectedTool === element ? "selected" : ''}
                                onClick={() => setSelectedTool(element)}
                            >
                                {element}
                            </button>
                        ))}
                        <button className={selectedTool === 'bond' ? 'selected' : ''} onClick={() => setSelectedTool('bond')}>Bond</button>
                        <button onClick={clearCanvas}>Clear</button>
                    </div>
                    <canvas ref={canvasRef} width="800" height="600" onClick={handleCanvasClick} className="drawing-area" />
                </div>
            </div>
        </div>
    );
}

export default App;
