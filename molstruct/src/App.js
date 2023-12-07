import React, { useState, useRef, useEffect } from 'react';
import './App.css';

function App() {
    const [selectedTool, setSelectedTool] = useState('C');
    const [atoms, setAtoms] = useState([]);
    const [bonds, setBonds] = useState([]);
    const [selectedAtoms, setSelectedAtoms] = useState([]);
    const canvasRef = useRef(null);
    const atomRadius = 8;

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

            /*
            // Debug circle
            context.beginPath();
            context.arc(x, y, atomRadius, 0, 2 * Math.PI);
            context.stroke();

            // Debug dot
            context.beginPath();
            context.arc(x, y, 2, 0, 2 * Math.PI);
            context.fillStyle = 'red';
            context.fill();
            */
        });

        // Draw bonds
        bonds.forEach(bond => {
            context.beginPath();
            context.moveTo(bond.start.x, bond.start.y);
            context.lineTo(bond.end.x, bond.end.y);
            context.stroke();
        });
    }, [atoms, bonds]);

    const handleCanvasClick = (event) => {
        const rect = canvasRef.current.getBoundingClientRect();
        const x = event.clientX - rect.left;
        const y = event.clientY - rect.top;

        const clickedAtom = atoms.find(atom =>
            Math.abs(atom.x - x) < 10 && Math.abs(atom.y - y) < 10);

        if (['C', 'N', 'O', 'H'].includes(selectedTool)) {
            setAtoms([...atoms, { type: selectedTool, x, y }]);
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
    }

    return (
        <div className="App">
            <div className="toolbar">
                <button onClick={() => setSelectedTool('C')}>C</button>
                <button onClick={() => setSelectedTool('N')}>N</button>
                <button onClick={() => setSelectedTool('O')}>O</button>
                <button onClick={() => setSelectedTool('H')}>H</button>
                <button onClick={() => setSelectedTool('bond')}>Bond</button>
            </div>
            <canvas ref={canvasRef} width="800" height="600" onClick={handleCanvasClick} className="drawing-area" />
        </div>
    );
}

export default App;
