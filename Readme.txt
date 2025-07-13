analogcircuitdsp
A powerful circuit simulation library built using Modified Nodal Analysis (MNA) for simulating electronic circuits in DC and transient (time-domain) regimes. It supports common linear components (resistors, capacitors, voltage/current sources) and essential non-linear components (diodes, MOSFETs, triodes, Op-Amps with saturation), and includes a functional delay line for audio effects.

Table of Contents
Introduction

Project Structure

Installation

Basic Usage

Creating a Circuit

Adding Components

Non-linear and Functional Components

DC Simulation

Transient Simulation

DelayLine Integration

API Reference

Extending the Library

Adding New Components

Creating New Complex Circuits

Contributing

License

1. Introduction
The analogcircuitdsp library provides Python-based tools to model and simulate electronic circuits accurately. It leverages the MNA method to construct the system of equations and scipy.optimize.fsolve to solve non-linear systems, making it suitable for simulating audio effects and other circuits with complex behaviors.

The primary goal of analogcircuitdsp is to offer a flexible and accurate framework for designers and developers to experiment with analog circuit designs in a digital environment, with a particular focus on faithfully simulating analog distortion and dynamic characteristics.

2. Project Structure
The analogcircuitdsp project is organized into the following main modules:

analogcircuitdsp/
├── components/                 # Definitions of individual electronic components
│   ├── __init__.py
│   ├── resistor.py             # Resistor class
│   ├── capacitor.py            # Capacitor class
│   ├── voltage_source.py       # VoltageSource class (can be time-dependent)
│   ├── current_source.py       # CurrentSource class (can be time-dependent)
│   ├── diode.py                # Diode class (non-linear)
│   ├── mosfet.py               # MOSFET class (non-linear)
│   ├── triode.py               # Triode/Vacuum Tube class (non-linear)
│   ├── op_amp.py               # OpAmp class (non-ideal, with saturation)
│   └── delay_line.py           # DelayLine class (functional, time-domain component)
├── circuits/                   # Definitions of pre-assembled / complex circuits
│   ├── __init__.py
│   ├── circuit.py              # Base class for defining any circuit
│   ├── passive_mixer.py        # Example: A simple passive resistive mixer
│   ├── gain_stage.py           # Example: An OpAmp-based gain stage
│   └── (your_custom_circuits).py
├── circuit_solver/             # MNA Solver modules
│   ├── __init__.py
│   └── mna_solver.py           # The core MNA solver logic
├── utils/                      # Utilities and helper functions
│   ├── __init__.py
│   └── helpers.py              # E.g., for numerical Jacobian (if not using fsolve's default)
├── main.py                     # Example usage or entry point
└── README.md                   # This file
3. Installation
To use the analogcircuitdsp library, you'll need Python 3.x installed. The primary dependencies are numpy and scipy.

Clone the repository:

Bash

git clone https://github.com/your_username/analogcircuitdsp.git
cd analogcircuitdsp
(Remember to replace your_username/analogcircuitdsp.git with the actual path to your repository.)

Install dependencies:

Bash

pip install numpy scipy
4. Basic Usage
The analogcircuitdsp library is designed to be modular. You typically start by creating an instance of a Circuit, then add components to it, define connections, and finally use the MnaSolver to run the simulation.

4.1. Creating a Circuit
All custom circuits you define should inherit from the analogcircuitdsp.circuits.circuit.Circuit class.

Python

# my_custom_circuit.py
from analogcircuitdsp.circuits.circuit import Circuit

class MySimpleCircuit(Circuit):
    def __init__(self, name="My Simple Circuit"):
        super().__init__(name)
        # Add circuit nodes. 'gnd' is the default 0V reference.
        self.add_node('input_node')
        self.add_node('output_node')
        self.add_node('gnd')
4.2. Adding Components
Import component classes from the analogcircuitdsp.components/ directory and add instances to your circuit using the add_component() method.

Python

# Continuing from my_custom_circuit.py
from analogcircuitdsp.components.resistor import Resistor
from analogcircuitdsp.components.voltage_source import VoltageSource

class MySimpleCircuit(Circuit):
    def __init__(self, name="My Simple Circuit"):
        super().__init__(name)
        self.add_node('input_node')
        self.add_node('output_node')
        self.add_node('gnd')
        
        # Add a DC voltage source at the input
        # Command: self.add_component(VoltageSource(name="component_name", node1="node_name1", node2="node_name2", voltage=value))
        self.add_component(VoltageSource(name="Vin", node1='input_node', node2='gnd', voltage=5.0))
        
        # Add a voltage divider using resistors
        # Command: self.add_component(Resistor(name="component_name", node1="node_name1", node2="node_name2", resistance=value))
        self.add_component(Resistor(name="R1", node1='input_node', node2='output_node', resistance=1000))
        self.add_component(Resistor(name="R2", node1='output_node', node2='gnd', resistance=2000))
4.3. Non-linear and Functional Components
Components like Diode, MOSFET, Triode, and OpAmp are automatically treated as non-linear by the MnaSolver. The DelayLine is a special functional component integrated during transient simulation.

Example: Adding an OpAmp to a GainStage

Python

# Example from analogcircuitdsp/circuits/gain_stage.py (or your custom circuit)
from analogcircuitdsp.components.op_amp import OpAmp
from analogcircuitdsp.components.resistor import Resistor
from analogcircuitdsp.circuits.circuit import Circuit

class GainStage(Circuit):
    def __init__(self, name="Gain_Stage", input_node='in', output_node='out', gain=2.0):
        super().__init__(name)
        
        # Public circuit nodes
        self.add_node(input_node)
        self.add_node(output_node)
        self.add_node('gnd')
        
        # Internal nodes for the OpAmp and feedback network
        self.add_node('opamp_inv')      # Inverting input of OpAmp
        self.add_node('opamp_noninv')   # Non-inverting input of OpAmp (connected to input_node)
        self.add_node('opamp_out_int')  # OpAmp's internal output node (before its Rout)

        # Calculate feedback resistor values for desired gain (e.g., non-inverting: Gain = 1 + Rf/Rg)
        Rf_val = (gain - 1) * 10e3 # Example calculation if Rg is 10k
        Rg_val = 10e3 
        
        # Instantiate and add the non-ideal OpAmp
        # Command: self.add_component(OpAmp(name="component_name", nodes={"in_inv": "node1", ...}, gain=value, ...))
        opamp_instance = OpAmp(name="U1", 
                               nodes={'in_inv': 'opamp_inv',
                                      'in_non_inv': input_node,
                                      'out': 'opamp_out_int'},
                               gain=1e5,             # High open-loop gain
                               input_resistance=1e12, # Very high differential input resistance
                               output_resistance=75.0, # Typical output resistance
                               saturation_voltage=13.0) # Saturation voltage for +/-15V rails
        self.add_component(opamp_instance)
        
        # IMPORTANT: Model the OpAmp's differential input resistance (Rin) explicitly as a Resistor.
        # The MNA solver treats this as a linear component.
        # Command: self.add_component(Resistor(name="component_name", node1="node1", node2="node2", resistance=value))
        self.add_component(Resistor(name="OpAmp_Rin_Model", 
                                    node1='opamp_inv', 
                                    node2=input_node, 
                                    resistance=opamp_instance.input_resistance))

        # Add feedback resistors for the OpAmp
        self.add_component(Resistor(name="Rf", node1='opamp_out_int', node2='opamp_inv', resistance=Rf_val))
        self.add_component(Resistor(name="Rg", node1='opamp_inv', node2='gnd', resistance=Rg_val))

        # Add an external load resistor from the OpAmp's output
        self.add_component(Resistor(name="R_load_out", node1='opamp_out_int', node2=output_node, resistance=10e3))

        # Map public nodes for easy access to input/output
        # Command: self.map_node("internal_node_name", "public_node_name")
        self.map_node(input_node, 'input')
        self.map_node(output_node, 'output')
Note on DelayLine Integration:
The analogcircuitdsp.components.delay_line.DelayLine is a functional component. It doesn't directly contribute to the MNA matrix equations. Instead, its behavior is integrated within the analogcircuitdsp.circuit_solver.mna_solver.MnaSolver.solve_transient loop. It takes a voltage from one circuit node as input and drives another node (typically via a VoltageSource) with its delayed output at each time step.

Python

# Example in a custom circuit using DelayLine
from analogcircuitdsp.components.delay_line import DelayLine
from analogcircuitdsp.components.voltage_source import VoltageSource
from analogcircuitdsp.components.resistor import Resistor
from analogcircuitdsp.circuits.circuit import Circuit 

class DelayEffectCircuit(Circuit):
    def __init__(self, max_delay_seconds=0.5, sample_rate=48000):
        super().__init__("Delay_Effect")
        self.add_node('in')
        self.add_node('gnd')
        self.add_node('delay_input_node')   # Node where signal enters the delay line
        self.add_node('delayed_output_node') # Node where delayed signal appears

        # Instantiate DelayLine as an attribute of the circuit instance
        # Command: self.delay_line = DelayLine(max_delay_seconds, sample_rate)
        self.delay_line = DelayLine(max_delay_seconds, sample_rate)
        
        # Define a VoltageSource to represent the output of the delay line.
        # The MnaSolver will update this source's voltage at each timestep using DelayLine.update().
        # Command: self.add_component(VoltageSource(name="component_name", node1="node1", node2="node2", voltage=initial_value))
        self.add_component(VoltageSource(name="V_delay_out_source", 
                                        node1='delayed_output_node', 
                                        node2='gnd', 
                                        voltage=0.0)) # Initial voltage
        
        # Connect circuit input to the delay_input_node (e.g., through a buffer resistor)
        self.add_component(Resistor(name="R_delay_in_buffer", node1='in', node2='delay_input_node', resistance=100.0))
        # Further circuit components would connect to 'delayed_output_node'
4.4. DC Simulation
To calculate the DC operating points (bias voltages) of the circuit:

Python

# main.py
from analogcircuitdsp.circuits.circuit import Circuit # Or your custom circuit class like MySimpleCircuit
from analogcircuitdsp.circuit_solver.mna_solver import MnaSolver
from analogcircuitdsp.components.resistor import Resistor
from analogcircuitdsp.components.voltage_source import VoltageSource
import numpy as np

# Example: Simple voltage divider circuit
class VoltageDivider(Circuit):
    def __init__(self):
        super().__init__("Voltage Divider")
        self.add_node('in')
        self.add_node('out')
        self.add_node('gnd')
        self.add_component(VoltageSource(name="Vin", node1='in', node2='gnd', voltage=10.0))
        self.add_component(Resistor(name="R1", node1='in', node2='out', resistance=1000))
        self.add_component(Resistor(name="R2", node1='out', node2='gnd', resistance=1000))

# 1. Create an instance of your circuit
my_circuit = VoltageDivider()

# 2. Create an MNA solver instance for the circuit
# Command: solver = MnaSolver(circuit_instance)
solver = MnaSolver(my_circuit)

# 3. Solve the circuit in DC
# Command: dc_solution_full, node_voltages, source_currents = solver.solve_dc(initial_guess=optional_array)
dc_solution_full, node_voltages, source_currents = solver.solve_dc(initial_guess=np.zeros(solver._total_equations))

if node_voltages is not None:
    print("DC Simulation completed:")
    # Retrieve node names and their corresponding voltages
    node_map = my_circuit.get_node_map()
    for node_name in node_map:
        node_idx = node_map[node_name]
        # Check if node_idx is valid (not -1 for 'gnd' if not explicitly added or connected)
        if node_idx != -1 and node_idx < len(node_voltages): 
            print(f"  Voltage at node '{node_name}': {node_voltages[node_idx]:.4f} V")
    # You can also print source currents if needed
    # print("\nVoltage Source Currents:")
    # for vs_name, vs_idx in solver._voltage_source_map.items():
    #     print(f"  Current through '{vs_name}': {source_currents[vs_idx]:.4f} A")
else:
    print("DC simulation did not converge.")
4.5. Transient Simulation
To simulate the circuit's evolution over time, for example, with audio signals:

Python

# main.py (continuation)
import matplotlib.pyplot as plt
import numpy as np
from analogcircuitdsp.circuits.circuit import Circuit
from analogcircuitdsp.components.voltage_source import VoltageSource
from analogcircuitdsp.components.resistor import Resistor
from analogcircuitdsp.circuit_solver.mna_solver import MnaSolver

# Example: Circuit with a time-dependent input
class ACInputCircuit(Circuit):
    def __init__(self, freq=1000, amplitude=1.0):
        super().__init__("AC_Test_Circuit")
        self.add_node('in')
        self.add_node('out')
        self.add_node('gnd')
        # Time-dependent voltage source (sine wave)
        self.add_component(VoltageSource(name="Vin", node1='in', node2='gnd', 
                                        voltage=lambda t: amplitude * np.sin(2 * np.pi * freq * t)))
        self.add_component(Resistor(name="R_load", node1='in', node2='out', resistance=100))
        self.add_component(Resistor(name="R_out", node1='out', node2='gnd', resistance=1000))

my_ac_circuit = ACInputCircuit(freq=440, amplitude=1.0)
solver = MnaSolver(my_ac_circuit)

# Transient simulation parameters
t_start = 0.0
t_end = 0.01    # Simulate for 10 ms
sample_rate = 48000 # Hz (This determines the simulation step: dt = 1/sample_rate)
dt = 1.0 / sample_rate

# Nodes to monitor (pass the public names if mapped, or internal names like 'input_node', 'output_node')
output_nodes_to_monitor = ['in', 'out'] 

# 1. Run the transient simulation
# Command: results, time_points = solver.solve_transient(start_time, end_time, time_step, output_nodes_list)
results, time_points = solver.solve_transient(t_start, t_end, dt, output_nodes=output_nodes_to_monitor)

if results is not None:
    print("Transient simulation completed.")
    # Plot results (requires matplotlib)
    plt.figure(figsize=(10, 6))
    plt.plot(results['time'], results['in'], label='Input Voltage (Node "in")')
    plt.plot(results['time'], results['out'], label='Output Voltage (Node "out")')
    plt.xlabel('Time (s)')
    plt.ylabel('Voltage (V)')
    plt.title(f'Transient Simulation of "{my_ac_circuit.name}"')
    plt.legend()
    plt.grid(True)
    plt.show()

    # For audio processing: the results in `results['out']` are your audio samples.
    # These can be saved to a WAV file or played back using appropriate libraries (e.g., `soundfile`, `pyaudio`).
else:
    print("Transient simulation failed to converge or encountered an error.")
4.6. DelayLine Integration
The DelayLine component is designed to be integrated during the solve_transient method. You instantiate it as an attribute of your Circuit class and update it within the simulation loop. The MnaSolver handles calling its update() method and applying its output to a VoltageSource in your circuit.

Key steps:

Instantiate DelayLine: Create self.delay_line = DelayLine(...) in your Circuit's __init__.

Define I/O Nodes: Create nodes like 'delay_input_node' and 'delayed_output_node' in your circuit.

Create Output VoltageSource: Add a VoltageSource to your circuit that will represent the DelayLine's output (e.g., VoltageSource(name="V_delay_out_source", node1='delayed_output_node', node2='gnd', voltage=0.0)). The MnaSolver will automatically find and update this source.

5. API Reference
For detailed API documentation, please refer to the source code docstrings. Below is a quick reference to key classes and methods, including their full import paths.

analogcircuitdsp.circuits.circuit.Circuit
Base class for defining any circuit.

__init__(self, name: str)

add_node(self, node_name: str)

add_component(self, component_instance): Accepts instances from analogcircuitdsp.components.

map_node(self, internal_node_name: str, public_node_name: str)

get_node_map(self) -> dict

get_num_nodes(self) -> int

get_num_voltage_sources(self) -> int

analogcircuitdsp.circuit_solver.mna_solver.MnaSolver
The core simulation engine.

__init__(self, circuit: Circuit)

solve_dc(self, initial_guess: np.ndarray = None) -> Tuple[np.ndarray, np.ndarray, np.ndarray]

solve_transient(self, t_start: float, t_end: float, dt: float, output_nodes: List[str] = None) -> Tuple[dict, np.ndarray]

analogcircuitdsp.components/ (Key Component Classes)
All component classes generally accept name (str), and node1, node2 (str), or a nodes dictionary (dict) for multi-terminal parts.

analogcircuitdsp.components.resistor.Resistor(name, node1, node2, resistance: float)

analogcircuitdsp.components.capacitor.Capacitor(name, node1, node2, capacitance: float)

analogcircuitdsp.components.voltage_source.VoltageSource(name, node1, node2, voltage: Union[float, Callable[[float], float]])

set_voltage(self, value: float): Updates fixed voltage.

get_voltage(self, time: float = 0) -> float: Gets current voltage (for time-dependent sources).

analogcircuitdsp.components.current_source.CurrentSource(name, node1, node2, current: Union[float, Callable[[float], float]])

set_current(self, value: float): Updates fixed current.

get_current(self, time: float = 0) -> float: Gets current (for time-dependent sources).

analogcircuitdsp.components.diode.Diode(name, node1, node2, is_current=1e-14, n_factor=1.0)

calculate_current(self, Vd: float) -> float

analogcircuitdsp.components.mosfet.MOSFET(name, nodes: dict, k_factor: float, threshold_voltage: float)

calculate_drain_current(self, Vgs: float, Vds: float) -> float

analogcircuitdsp.components.triode.Triode(name, nodes: dict, mu: float, kp: float, kv: float)

calculate_plate_current(self, Vgk: float, Vpk: float) -> float

analogcircuitdsp.components.op_amp.OpAmp(name, nodes: dict, gain=1e5, input_resistance=1e6, output_resistance=50.0, saturation_voltage=15.0)

calculate_output_voltage(self, v_plus: float, v_minus: float) -> float

get_input_resistance() -> float

get_output_resistance() -> float

analogcircuitdsp.components.delay_line.DelayLine(max_delay_seconds: float, sample_rate: float)

set_delay_time(self, delay_time_seconds: float)

update(self, input_sample: float) -> float

reset(self)

analogcircuitdsp.circuits/ (Pre-defined Circuit Examples)
These modules contain ready-to-use circuit implementations.

analogcircuitdsp.circuits.passive_mixer.PassiveMixer(name: str, num_inputs: int, mix_res_val: float)

analogcircuitdsp.circuits.gain_stage.GainStage(name: str, input_node: str, output_node: str, gain: float)

6. Extending the Library
analogcircuitdsp is designed for easy extension.

6.1. Adding New Components
To add a new type of electronic component:

Create a new Python file in the analogcircuitdsp/components/ directory (e.g., inductor.py).

Define a class for the component, including its name, connecting nodes (node1, node2, or a nodes dictionary for multi-terminal components), and any relevant parameters (e.g., inductance, specific SPICE parameters).

Implement any necessary methods required by the MnaSolver (e.g., methods to calculate currents, voltages, or contributions to the MNA matrices/Jacobian, if it's non-linear).

Update analogcircuitdsp/circuit_solver/mna_solver.py:

Add the import for your new component.

Modify the _build_linear_system method if it's a linear component (e.g., inductor in transient analysis).

Modify the _system_equations method if it's a non-linear component (it will need to contribute to the F vector).

6.2. Creating New Complex Circuits
To assemble more complex circuits from existing components:

Create a new Python file in the analogcircuitdsp/circuits/ directory (e.g., fuzz_pedal.py).

Define a class that inherits from analogcircuitdsp.circuits.circuit.Circuit.

Inside your class's __init__ method:

Define all necessary internal nodes using self.add_node().

Instantiate and add all individual components (Resistors, Capacitors, Diodes, OpAmps, etc.) using self.add_component(), connecting them to your defined nodes.

If your circuit has "public" input/output nodes, use self.map_node() to make them easily accessible from outside the circuit definition.

If your circuit includes functional components like DelayLine, add them as attributes of your Circuit class instance, not as add_component() calls.

7. Contributing
Contributions to analogcircuitdsp are highly encouraged! If you find a bug, have a suggestion for a new feature, or want to contribute code, feel free to open an issue or submit a pull request.

8. License
This project is licensed under the [Your License Here, e.g., MIT License]. See the LICENSE file for more details.







