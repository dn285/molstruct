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
    const bondLength = 50;
    const bondSpace = 2; // The space between double bonds

    const [structuralData, setStructuralData] = useState({
        molecular: '',
        smiles: '',
        stdinchi: '',
        stdinchikey: '',
        iupac: ''
    });

    const clearCanvas = () => {
        setAtoms([]);
        setBonds([]);
    };

    useEffect(() => {

        // Bond drawing functions
        function drawSingleBond(context, startX, startY, endX, endY) {
            // Calculate direction vector
            const dx = endX - startX;
            const dy = endY - startY;

            const dist = Math.sqrt(dx * dx + dy * dy);

            const vx = dx / dist;
            const vy = dy / dist;

            const offsetX = vx * atomRadius;
            const offsetY = vy * atomRadius;

            context.beginPath();
            context.moveTo(startX + offsetX, startY + offsetY);
            context.lineTo(endX - offsetX, endY - offsetY);
            context.stroke();
        }

        function drawDoubleBond(context, startX, startY, endX, endY, bondSpace) {
            // Calculate direction vector
            const dx = endX - startX;
            const dy = endY - startY;

            // Calculate perpendicular vector
            const dist = Math.sqrt(dx * dx + dy * dy);
            const ux = -dy / dist;
            const uy = dx / dist;

            // Calculate offset points
            const sx1 = startX + ux * bondSpace;
            const sy1 = startY + uy * bondSpace;
            const ex1 = endX + ux * bondSpace;
            const ey1 = endY + uy * bondSpace;

            const sx2 = startX - ux * bondSpace;
            const sy2 = startY - uy * bondSpace;
            const ex2 = endX - ux * bondSpace;
            const ey2 = endY - uy * bondSpace;

            drawSingleBond(context, sx1, sy1, ex1, ey1);
            drawSingleBond(context, sx2, sy2, ex2, ey2);
        }

        function drawTripleBond(context, startX, startY, endX, endY) {
            drawSingleBond(context, startX, startY, endX, endY);
            drawDoubleBond(context, startX, startY, endX, endY, 2 * bondSpace);
        }

        // Interface with backend
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

        // Drawing the canvas

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

            // Render implicit hydrogens 
            const hydrogens = calculateHydrogens(atom, bonds);
            const hydrogenPositions = calculateHydrogenPositions(atom, hydrogens, bondLength)
            hydrogenPositions.forEach(pos => {
                context.font = "16px Arial";
                context.textAlign = "center";
                context.textBaseline = "middle";
                context.fillText('H', pos.x, pos.y);

                drawSingleBond(context, x, y, pos.x, pos.y);
            });
        });

        // Draw bonds
        bonds.forEach(bond => {
            const startAtom = bond.start;
            const endAtom = bond.end;

            if (bond.multiplicity === 1) {
                drawSingleBond(context, startAtom.x, startAtom.y, endAtom.x, endAtom.y);
            }
            else if (bond.multiplicity === 2) {
                drawDoubleBond(context, startAtom.x, startAtom.y, endAtom.x, endAtom.y, bondSpace);
            }
            else if (bond.multiplicity === 3) {
                drawTripleBond(context, startAtom.x, startAtom.y, endAtom.x, endAtom.y);
            }
        });

        updateStructuralData();
    }, [atoms, bonds]);

    // Implicit hydrogen calculations

    const calculateHydrogens = (atom, bonds) => {
        const valence = {
            'B': 3, 'C': 4, 'N': 3, 'O': 2, 'P': 3
        };
        const bondCount = bonds
            .filter(bond => bond.start.id === atom.id || bond.end.id === atom.id)
            .reduce((sum, bond) => sum + bond.multiplicity, 0);

        const hydrogensNeeded = valence[atom.type] - bondCount;
        return hydrogensNeeded > 0 ? hydrogensNeeded : 0;
    };

    function calculateHydrogenPositions(atom, hydrogenCount, bondLength) {
        const centerX = atom.x;
        const centerY = atom.y;

        // Calculate existing bond angles around this atom
        const existingAngles = bonds.filter(bond => bond.start.id === atom.id || bond.end.id === atom.id)
            .map(bond => {
                const bondedAtom = bond.start.id === atom.id ? bond.end : bond.start;
                return Math.atan2(bondedAtom.y - centerY, bondedAtom.x - centerX);
            });

        // Find largest gap
        let largestGap = 0;
        let largestStart = 0;

        if (existingAngles.length <= 1) {
            largestGap = 2 * Math.PI;
            largestStart = existingAngles.length === 1 ? existingAngles[0] : 0;
        }
        else {
            for (let i = 0; i < existingAngles.length; i++) {
                const currentAngle = existingAngles[i];
                const nextAngle = existingAngles[(i + 1) % existingAngles.length];
                const gap = (nextAngle - currentAngle + 2 * Math.PI) % (2 * Math.PI);

                if (gap > largestGap) {
                    largestGap = gap;
                    largestStart = currentAngle;
                }
            }
        }

        // Evenly place hydrogens into largest gap
        const positions = [];
        const angleInterval = largestGap / (hydrogenCount + (existingAngles.length > 0));
        for (let i = 1; i <= hydrogenCount; i++) {
            const angle = largestStart + angleInterval * i;
            const x = centerX + bondLength * Math.cos(angle);
            const y = centerY + bondLength * Math.sin(angle);
            positions.push({ x, y });
        }

        return positions;
    }

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
        else {
            if (clickedAtom) {
                if (selectedAtoms.length < 1) {
                    setSelectedAtoms([clickedAtom]);
                }
                else {
                    if (selectedAtoms[0].id !== clickedAtom.id) {
                        const multiplicity = Number(selectedTool.slice(-1));
                        const newBond = createBond(selectedAtoms[0], clickedAtom, multiplicity);
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

    const createBond = (startAtom, endAtom, multiplicity) => {
        return {
            start: {
                id: startAtom.id,
                x: startAtom.x,
                y: startAtom.y
            },
            end: {
                id: endAtom.id,
                x: endAtom.x,
                y: endAtom.y
            },
            multiplicity: multiplicity
        };
    };

    // The actual page

    return (
        <div className="App">
            <div className="top-container">
                <div className="sidebar">
                    <h1>Welcome to molstruct!</h1>
                    <p>Draw your molecule onto the canvas on the right and see
                        the corresponding structural data and IUPAC name below.
                        Click to place atoms, and click two atoms to draw a
                        bond between them. Note that IUPAC names are not
                        guaranteed to be correct.</p>
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
                        <button className={selectedTool === 'bond1' ? 'selected' : ''} onClick={() => setSelectedTool('bond1')}>Bond1</button>
                        <button className={selectedTool === 'bond2' ? 'selected' : ''} onClick={() => setSelectedTool('bond2')}>Bond2</button>
                        <button className={selectedTool === 'bond3' ? 'selected' : ''} onClick={() => setSelectedTool('bond3')}>Bond3</button>
                        <button onClick={clearCanvas}>Clear</button>
                    </div>
                    <canvas ref={canvasRef} width="800" height="600" onClick={handleCanvasClick} className="drawing-area" />
                </div>
            </div>
        </div>
    );
}

export default App;
