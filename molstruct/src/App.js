import React, { useState, useRef, useEffect } from 'react';
import './App.css';

function App() {
    const [selectedTool, setSelectedTool] = useState('C');
    const [atoms, setAtoms] = useState([]);
    const [bonds, setBonds] = useState([]);
    const [selectedAtoms, setSelectedAtoms] = useState([]);
    const canvasRef = useRef(null);

    useEffect(() => {
        const canvas = canvasRef.current;
        const context = canvas.getContext('2d');

        // Clear canvas
        context.clearRect(0, 0, canvas.width, canvas.height);

        // Draw atoms
        atoms.forEach(atom => {
            context.font = "20px Arial";
            context.fillText(atom.type, atom.x, atom.y);
        })

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
                    setBonds([...bonds, { start: selectedAtoms[0], end: clickedAtom }])
                    setSelectedAtoms([]);
                }
            }
        }
    };

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
