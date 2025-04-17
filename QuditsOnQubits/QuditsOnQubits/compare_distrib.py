from qiskit.quantum_info import Statevector
from qiskit_aer import AerSimulator
from qiskit_ibm_runtime import SamplerV2
from qiskit_ibm_runtime.fake_provider import FakeTorino
from qiskit_ibm_runtime.fake_provider import FakeBrisbane
import matplotlib.pyplot as plt
import numpy as np

from QuditsOnQubits import torino, brisbane

def compare_distrb(circ):
    torinocirc = torino(circ)
    briscirc = brisbane(circ)

    circ_temp = circ.copy()
    circ_temp.remove_final_measurements()
    d = Statevector.from_instruction(circ_temp).probabilities_dict()

    aerBris = AerSimulator.from_backend(FakeBrisbane())
    samplerBris = SamplerV2(aerBris)
    jobBris = samplerBris.run([briscirc])
    pub_resultBris = jobBris.result()[0]
    countsBris = pub_resultBris.data.meas.get_counts()

    aerTorino = AerSimulator.from_backend(FakeTorino())
    samplerTorino = SamplerV2(aerTorino)
    jobTorino = samplerTorino.run([torinocirc])
    pub_resultTorino = jobTorino.result()[0]
    countsTorino = pub_resultTorino.data.meas.get_counts()

    keys = np.array(list(d.keys()))
    vals = np.array(list(d.values())) * 1000

    keysBris = list(countsBris.keys())
    valsBris = list(countsBris.values())

    keysTorino = list(countsTorino.keys())
    valsTorino = list(countsTorino.values())

    plt.bar(keys, vals, color="blue", label='Pure')
    plt.bar(keysBris, valsBris, color="red", label='Brisbane', alpha=0.7)
    plt.bar(keysTorino, valsTorino, color="green", label='Torino', alpha=0.7)

    plt.xticks([])
    plt.legend()
    plt.show()