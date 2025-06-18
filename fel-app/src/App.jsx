import React from 'react';
import './App.css';
import Dropdown from './components/Dropdown/Dropdown';
function App()
{
    return (
        <>
        <div className="layout">
          <div className="sidebar">
            <h2>Menu</h2>
            <Dropdown buttonText="test" contentText = <p>123</p>/>
          </div>
          <div className="main-content">
              <h1>FEL simulation</h1>
              <p>select option on left</p>
          </div>
          <div className="linegraph ">
            <h1>graph here</h1>
          </div>
          <div className="twiss">
            <h10>Twiss options</h10>
          </div>
        </div>
        </>
    );
}

export default App;
