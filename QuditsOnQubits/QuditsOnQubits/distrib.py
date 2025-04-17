from qiskit.quantum_info import Statevector
from qiskit.visualization import plot_histogram
from qiskit_aer import AerSimulator
from qiskit_ibm_runtime import SamplerV2
from qiskit_ibm_runtime.fake_provider import FakeTorino
from qiskit_ibm_runtime.fake_provider import FakeBrisbane
import matplotlib.pyplot as plt
from qiskit.visualization import plot_histogram
import numpy as np

def distrib(circ, plot=True):

    if circ.name != "brisbane" and circ.name != "torino":
        circ_temp = circ.copy()
        circ_temp.remove_final_measurements()
        d = Statevector.from_instruction(circ_temp).probabilities_dict()

        if plot:
            keys = np.array(list(d.keys()))
            vals = np.array(list(d.values())) * 1000
            plt.bar(keys, vals)
            plt.xticks([])
            plt.show()
            return d
        else:
            return d

    else:
        if circ.name == "brisbane":
            backend = FakeBrisbane()
        elif circ.name == "torino":
            backend = FakeTorino()
        aer = AerSimulator.from_backend(backend)
        sampler = SamplerV2(aer)
        job = sampler.run([circ])
        pub_result = job.result()[0]
        counts = pub_result.data.meas.get_counts()
        if plot:
            plot_histogram(counts)
            return counts
        else:
            return counts