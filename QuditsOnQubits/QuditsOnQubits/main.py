import numpy as np
from qiskit.circuit import QuantumCircuit
from qiskit.quantum_info import Statevector, Operator
from qiskit.synthesis import TwoQubitWeylDecomposition
from qiskit.circuit.library import UnitaryGate
from qiskit_ibm_runtime.fake_provider import FakeBrisbane
from qiskit_ibm_runtime.fake_provider import FakeTorino
from qiskit.transpiler.preset_passmanagers import generate_preset_pass_manager
from qiskit_ibm_runtime import SamplerV2
from qiskit_aer import AerSimulator
from qiskit.visualization import plot_histogram
from qiskit.quantum_info import hellinger_fidelity
import matplotlib.pyplot as plt

from QuditsOnQubits.ghz_creator import create_basic_ghz
from QuditsOnQubits.ghz_creator import create_graph_ghz

class QuditsOnQubits:
    def __init__(self, graph=None):
        if graph is None:
            self.circuit = create_basic_ghz()

        else:
            self.circuit = create_graph_ghz(graph)

        self.backend = None
        self.transpiled = None
        self.statBasic = None
        self.statTorino = None
        self.statBrisbane = None
        self.name = None
    def brisbane(self):
        backendBrisbane = FakeBrisbane()
        self.backend = backendBrisbane
        self.name = "Brisbane"
        targetBrisbane = backendBrisbane.target

        pass_managerBrisbane = generate_preset_pass_manager(
            optimization_level=3, target=targetBrisbane
        )

        transpiledBrisbane = pass_managerBrisbane.run(self.circuit)

        self.transpiled = transpiledBrisbane

        return transpiledBrisbane


    def torino(self):

        backendTorino = FakeTorino()
        self.backend = backendTorino
        self.name = "Torino"

        targetTorino = backendTorino.target

        pass_managerTorino = generate_preset_pass_manager(
            optimization_level=3, target=targetTorino
        )

        transpiledTorino = pass_managerTorino.run(self.circuit)

        self.transpiled = transpiledTorino

        return transpiledTorino

    def distrb(self):
        if self.backend is None:
            d = Statevector.from_instruction(self.circuit).probabilities_dict()
            self.statBasic = d
            plot_histogram(d, figsize=(30,10))
        else:
            aer = AerSimulator.from_backend(self.backend)
            sampler = SamplerV2(aer)
            job = sampler.run([self.transpiled])
            pub_result = job.result()[0]
            counts = pub_result.data.meas.get_counts()
            if self.name == "Brisbane":
                self.statBrisbane = counts
            elif self.name == "Torino":
                self.statTorino = counts
            plot_histogram(counts, figsize=(30, 10)).show()
    def fidelity(self, backend_name):
        if backend_name == "Brisbane":
            fid = hellinger_fidelity(self.statBasic, self.statBrisbane)
        elif backend_name == "Torino":
            fid = hellinger_fidelity(self.statBasic, self.statTorino)

        return fid