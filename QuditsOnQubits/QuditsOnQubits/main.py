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
        self.name = None
        tempCirc = self.circuit.copy()
        tempCirc.remove_final_measurements()
        self.statBasic = Statevector.from_instruction(tempCirc).probabilities_dict()

        self.statTorino = None
        self.statBrisbane = None

    def brisbane(self):
        backend_brisbane = FakeBrisbane()
        self.backend = backend_brisbane
        self.name = "Brisbane"


        target_brisbane = backend_brisbane.target
        pass_manager_brisbane = generate_preset_pass_manager(
            optimization_level=3, target=target_brisbane
        )
        transpiled_brisbane = pass_manager_brisbane.run(self.circuit)
        self.transpiled = transpiled_brisbane


        aer = AerSimulator.from_backend(backend_brisbane)
        sampler = SamplerV2(aer)
        job = sampler.run([transpiled_brisbane])
        pub_result = job.result()[0]
        counts = pub_result.data.meas.get_counts()


        self.statBrisbane = counts

        return transpiled_brisbane


    def torino(self):

        backend_torino = FakeTorino()
        self.backend = backend_torino
        self.name = "Torino"

        target_torino = backend_torino.target
        pass_manager_torino = generate_preset_pass_manager(
            optimization_level=3, target=target_torino
        )
        transpiled_torino = pass_manager_torino.run(self.circuit)
        self.transpiled = transpiled_torino


        aer = AerSimulator.from_backend(backend_torino)
        sampler = SamplerV2(aer)
        job = sampler.run([transpiled_torino])
        pub_result = job.result()[0]
        counts = pub_result.data.meas.get_counts()

        self.statTorino = counts

        return transpiled_torino

    def distrb(self):
        if self.backend is None:
            circ_temp = self.circuit.copy()
            circ_temp.remove_final_measurements()
            d = Statevector.from_instruction(circ_temp).probabilities_dict()
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
            plot_histogram(counts, figsize=(30, 10))
    def fidelity(self, backend_name):
        if backend_name == "Brisbane":
            if self.statBrisbane is None:
                raise ValueError("First use .brisbane() to count fidelity.")
            fid = hellinger_fidelity(self.statBasic, self.statBrisbane)

        elif backend_name == "Torino":
            if self.statTorino is None:
                raise ValueError("First use .torino() to count fidelity.")
            fid = hellinger_fidelity(self.statBasic, self.statTorino)

        else:
            raise ValueError("Unknown backend. Use \"Brisbane\" or \"Torino\".")

        return fid