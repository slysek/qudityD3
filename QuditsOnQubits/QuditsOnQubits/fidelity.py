from qiskit.quantum_info import hellinger_fidelity
from QuditsOnQubits import distrib
def fidelity(base, compare):

    counts_base = distrib(base, plot=False)
    counts_compare = distrib(compare, plot=False)

    fid = hellinger_fidelity(counts_base, counts_compare)
    return fid