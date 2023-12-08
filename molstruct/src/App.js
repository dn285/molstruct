import React, { useState, useRef, useEffect } from 'react';
import './App.css';

function App() {
    const [selectedTool, setSelectedTool] = useState('C');
    const [atoms, setAtoms] = useState([]);
    const [bonds, setBonds] = useState([]);
    const [selectedAtoms, setSelectedAtoms] = useState([]);
    const canvasRef = useRef(null);
    const atomRadius = 8;

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
    }, [atoms, bonds]);

    // Handlers

    const handleCanvasClick = (event) => {
        const rect = canvasRef.current.getBoundingClientRect();
        const x = event.clientX - rect.left;
        const y = event.clientY - rect.top;

        const clickedAtom = atoms.find(atom =>
            Math.abs(atom.x - x) < 10 && Math.abs(atom.y - y) < 10);

        if (['C', 'N', 'O', 'H'].includes(selectedTool)) {
            addAtom(selectedTool, x, y);
            setSelectedAtoms([]);
        }
        else if (selectedTool == 'bond') {
            if (clickedAtom) {
                if (selectedAtoms.length < 1) {
                    setSelectedAtoms([clickedAtom]);
                }
                else {
                    const newBond = createBond(selectedAtoms[0], clickedAtom, atomRadius);
                    setBonds([...bonds, newBond])
                    setSelectedAtoms([]);
                }
            }
        }
    };

    const addAtom = (atomType, x, y) => {
        const newAtom = { type: atomType, x: x, y: y };
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
                x: startAtom.x + offsetX,
                y: startAtom.y + offsetY
            },
            end: {
                x: endAtom.x - offsetX,
                y: endAtom.y - offsetY
            }
        };
    };

    // Computing formulae

    const calculateMolecularFormula = (atoms) => {
        const atomCounts = atoms.reduce((acc, atom) => {
            acc[atom.type] = (acc[atom.type] || 0) + 1;
            return acc;
        }, {});

        let formula = '';
        for (const [atom, count] of Object.entries(atomCounts)) {
            formula += `${atom}${count > 1 ? count : ''}`;
        }
        return formula;
    };

    // The actual page

    return (
        <div className="App">
            <div className="top-container">
                <div className="sidebar">
                    <h1>Welcome to molstruct!</h1>
                    <p>Draw your molecule onto the canvas on the right and see the corresponding IUPAC name below.</p>
                    <div className="output-names">
                        <h2>Outputs</h2>
                        <div className="name-group">
                            <label>Molecular Formula</label>
                            <input type="text" readOnly value={calculateMolecularFormula(atoms)} />
                        </div>
                        <div className="name-group">
                            <label>IUPAC Name</label>
                            <input type="text" readOnly value="" />
                        </div>
                    </div>
                </div>
                <div className="canvas-container">
                    <div className="toolbar">
                        <button className={selectedTool === 'C' ? 'selected' : ''} onClick={() => setSelectedTool('C')}>C</button>
                        <button className={selectedTool === 'N' ? 'selected' : ''} onClick={() => setSelectedTool('N')}>N</button>
                        <button className={selectedTool === 'O' ? 'selected' : ''} onClick={() => setSelectedTool('O')}>O</button>
                        <button className={selectedTool === 'H' ? 'selected' : ''} onClick={() => setSelectedTool('H')}>H</button>
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
