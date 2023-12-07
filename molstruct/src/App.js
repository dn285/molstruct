import React, { useState, useRef, useEffect } from 'react';
import './App.css';

function App() {
    const [selectedTool, setSelectedTool] = useState('dot');
    const [dots, setDots] = useState([]);
    const [lines, setLines] = useState([]);
    const [selectedDots, setSelectedDots] = useState([]);
    const canvasRef = useRef(null);

    useEffect(() => {
        const canvas = canvasRef.current;
        const context = canvas.getContext('2d');

        // Clear canvas
        context.clearRect(0, 0, canvas.width, canvas.height);

        // Draw dots
        dots.forEach(dot => {
            context.beginPath();
            context.arc(dot.x, dot.y, 5, 0, 2 * Math.PI);
            context.fill();
        });

        // Draw lines
        lines.forEach(line => {
            context.beginPath();
            context.moveTo(line.start.x, line.start.y);
            context.lineTo(line.end.x, line.end.y);
            context.stroke();
        });
    }, [dots, lines]);

    const handleCanvasClick = (event) => {
        const rect = canvasRef.current.getBoundingClientRect();
        const x = event.clientX - rect.left;
        const y = event.clientY - rect.top;
        const clickedDot = dots.find(dot => Math.hypot(dot.x - x, dot.y - y) < 5);

        if (selectedTool === 'dot') {
            setDots([...dots, { x, y }]);
        }
        else if (selectedTool == 'line') {
            if (clickedDot) {
                if (selectedDots.length < 1) {
                    setSelectedDots([clickedDot]);
                }
                else {
                    setLines([...lines, { start: selectedDots[0], end: clickedDot }])
                    setSelectedDots([]);
                }
            }
        }
    };

    return (
        <div className="App">
            <div className="toolbar">
                <button onClick={() => setSelectedTool('dot')}>Dot Tool</button>
                <button onClick={() => setSelectedTool('line')}>Line Tool</button>
            </div>
            <canvas ref={canvasRef} width="800" height="600" onClick={handleCanvasClick} className="drawing-area" />
        </div>
    );
}

export default App;
