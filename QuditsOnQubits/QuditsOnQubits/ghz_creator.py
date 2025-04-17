import numpy as np
from qiskit.circuit import QuantumCircuit
from qiskit.quantum_info import Statevector, Operator
from qiskit.circuit.library import UnitaryGate
from igraph import Graph
from qiskit import qpy
def create_ghz(raw=False):

    with open('GHZ.qpy', 'rb') as fd:
        ghz = qpy.load(fd)[0]

    fgate = get_f_gate()

    if raw == True:
        return ghz, Statevector(ghz)
    else:
        ghz_circ = ghz.copy()
        ghz_circ.remove_final_measurements()
        cartanH = Operator(fgate).conjugate().transpose()
        changeToDefault = cartanH.tensor(cartanH).tensor(Operator(np.identity(4))).to_matrix()
        optiGhzState = Statevector(ghz_circ)
        optiGhzMatrix = optiGhzState.data.reshape((64,1))

        default = np.dot(changeToDefault, optiGhzMatrix)

        return ghz, Statevector(default)

def get_f_gate():
    with open('Fgate.qpy', 'rb') as fd:
        new_qc = qpy.load(fd)[0]
    return new_qc

def create_graph_ghz(graph):
    with open('CZ.qpy', 'rb') as fd:
        CZ = qpy.load(fd)[0]

    with open('Fgate.qpy', 'rb') as fd:
        cartanCircuit = qpy.load(fd)[0]
    def createGraphCirc(G, cartanCirc):
        n = G.vcount()  # liczba wierzchołków w igraph
        qubitList = [[2 * i, 2 * i + 1] for i in range(n)]

        edgeList = []

        # W igraph kolejność krawędzi jest zachowana domyślnie
        for edge in G.get_edgelist():
            u, v = edge[0], edge[1]
            if u != v:
                edgeList.append([
                    qubitList[u][0],
                    qubitList[u][1],
                    qubitList[v][0],
                    qubitList[v][1]
                ])

        graphCirc = QuantumCircuit(2 * n)

        for pair in qubitList:
            graphCirc.append(cartanCirc, pair)

        # Dodajemy bramki CZ w kolejności zgodnej z kolejnością dodania krawędzi w igraph
        for edge in edgeList:
            graphCirc.append(CZ, edge)

        return graphCirc


    basisTransformationMatrix = np.identity(2)
    uniGate = UnitaryGate(basisTransformationMatrix)

    final_circ = createGraphCirc(graph, cartanCircuit)

    return final_circ

