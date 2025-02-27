import numpy as np
from qiskit.circuit import QuantumCircuit
from qiskit.quantum_info import Statevector, Operator
from qiskit.synthesis import TwoQubitWeylDecomposition
from qiskit.circuit.library import UnitaryGate
import networkx as nx
def create_basic_ghz():

    #Creating F component of graph
    def cartanCircuit(VCxGate, transGate):
        circ = QuantumCircuit(2)
        phi = np.arctan(1 + np.sqrt(3))


        #qubit 0 after VCx
        circ.append(transGate, [0])
        circ.append(transGate, [1])
        circ.z(0)
        circ.x(0)
        circ.ry(phi, 0)
        circ.x(0)
        circ.p(15 * np.pi / 16, 0)
        circ.x(0)
        circ.p(3 * np.pi / 16, 0)
        circ.x(0)

        #qubit 1 after VCx
        circ.z(1)
        circ.x(1)
        circ.ry(phi, 1)
        circ.x(1)
        circ.z(1)
        circ.x(1)
        circ.p(13 * np.pi / 16, 1)
        circ.x(1)
        circ.p(9 * np.pi / 16, 1)

        circ.barrier(0)
        circ.barrier(1)

        circ.append(VCxGate.decompose(), [0, 1])


        #qubit 0
        circ.p(15 * np.pi / 16, 0)
        circ.x(0)
        circ.p(3 * np.pi / 16, 0)
        circ.x(0)
        circ.z(0)
        circ.x(0)
        circ.ry(phi, 0)
        circ.x(0)

        #qubit 1
        circ.p(15 * np.pi / 16, 1)
        circ.x(1)
        circ.p(3 * np.pi / 16, 1)
        circ.x(1)
        circ.z(1)
        circ.x(1)
        circ.ry(phi, 1)
        circ.x(1)

        circ.append(transGate.conjugate(), [0])
        circ.append(transGate.conjugate(), [1])

        return circ

    pi = np.pi
    sqrt = np.sqrt
    log = np.log

    common_log_term = log(-1 / 6 * (-1) ** (1 / 4) * (sqrt(24 - 6 * sqrt(3)) + 1j * (3 + sqrt(3))))
    common_factor = 3j + 1j * sqrt(3) + sqrt(6 * (4 - sqrt(3)))
    denominator = 15 * pi - 4 * 1j * common_log_term

    #Matrix from cartan decomposition
    VCx = np.array([
        [
            -(((-1) ** (5 / 8) * common_factor * (30 * pi - 8 * 1j * common_log_term)) / (24 * denominator)) +
            (3 * (-1) ** (5 / 8) * (-30 * pi + 8 * 1j * common_log_term)) / (2 * common_factor * denominator),
            0,
            0,
            (3 * (-1) ** (5 / 8) * (-30 * pi + 8 * 1j * common_log_term)) / (2 * common_factor * denominator) -
            ((-1) ** (5 / 8) * common_factor * (-30 * pi + 8 * 1j * common_log_term)) / (24 * denominator)
        ],
        [
            0,
            -(1 / 2) * (-1) ** (1 / 8) - (1 / 2) * (-1) ** (5 / 8),
            -(1 / 2) * (-1) ** (1 / 8) + (1 / 2) * (-1) ** (5 / 8),
            0
        ],
        [
            0,
            -(1 / 2) * (-1) ** (1 / 8) + (1 / 2) * (-1) ** (5 / 8),
            -(1 / 2) * (-1) ** (1 / 8) - (1 / 2) * (-1) ** (5 / 8),
            0
        ],
        [
            (3 * (-1) ** (5 / 8) * (-30 * pi + 8 * 1j * common_log_term)) / (2 * common_factor * denominator) -
            ((-1) ** (5 / 8) * common_factor * (-30 * pi + 8 * 1j * common_log_term)) / (24 * denominator),
            0,
            0,
            -(((-1) ** (5 / 8) * common_factor * (30 * pi - 8 * 1j * common_log_term)) / (24 * denominator)) +
            (3 * (-1) ** (5 / 8) * (-30 * pi + 8 * 1j * common_log_term)) / (2 * common_factor * denominator)
        ]
    ], dtype=complex)

    #Creating circuit from matrix using Weyl decomposition
    weylDecomCirc = TwoQubitWeylDecomposition(VCx)
    VCxGate = weylDecomCirc.circuit()

    basisTransformationMatrix = np.identity(2)
    uniGate = UnitaryGate(basisTransformationMatrix)
    cartanCirc = cartanCircuit(VCxGate, uniGate)

    #CZ gate circuit
    def CZqutrit(transGate):
        circ = QuantumCircuit(4, name="CZ")

        circ.append(transGate, [0])
        circ.append(transGate, [1])
        circ.append(transGate, [2])
        circ.append(transGate, [3])
        omega = 2*np.pi/3

        circ.cp(omega, 0, 2)
        circ.cp(omega, 1, 3)
        circ.cp(-1 * omega, 0, 3)
        circ.cp(-1 * omega, 1, 2)

        circ.append(transGate.conjugate().transpose(), [0])
        circ.append(transGate.conjugate().transpose(), [1])
        circ.append(transGate.conjugate().transpose(), [2])
        circ.append(transGate.conjugate().transpose(), [3])

        return circ

    czqutrit = CZqutrit(uniGate)

    #Making ghz circuit using caratinCirc (F operator) and czqutrit (CZ operator)
    ghzcirc = QuantumCircuit(6)
    ghzcirc.append(cartanCirc, [0, 1])
    ghzcirc.append(cartanCirc, [2, 3])
    ghzcirc.append(cartanCirc, [4, 5])
    ghzcirc.append(czqutrit, [0, 1, 2, 3])
    ghzcirc.append(czqutrit, [0, 1, 4, 5])

    #Optimizing F gate
    ghzState = Statevector(ghzcirc)
    ghzMatrix = ghzState.data.reshape((64,1))
    testWeyly = TwoQubitWeylDecomposition(Operator(cartanCirc).to_matrix())
    optiF = testWeyly.circuit()
    optiF.name = "+"

    #Optimizing ghz circuit
    optiGHZ = QuantumCircuit(6)
    optiGHZ.append(optiF, [0, 1])
    optiGHZ.append(optiF, [2, 3])
    optiGHZ.append(optiF, [4, 5])
    optiGHZ.append(czqutrit, [0, 1, 2, 3])
    optiGHZ.append(czqutrit, [0, 1, 4, 5])
    optiGHZ.measure_all()

    #%%

    ghz_circ = optiGHZ.copy()
    ghz_circ.remove_final_measurements()
    cartanH = Operator(cartanCirc).conjugate().transpose()
    changeToDefault = cartanH.tensor(cartanH).tensor(Operator(np.identity(4))).to_matrix()
    optiGhzState = Statevector(ghz_circ)
    optiGhzMatrix = optiGhzState.data.reshape((64,1))

    default = np.dot(changeToDefault, optiGhzMatrix)

    return [optiGHZ, default]

def create_graph_ghz(graph):
    def cartanCircuit(VCxGate, transGate):
        circ = QuantumCircuit(2, name='+')
        phi = np.arctan(1 + np.sqrt(3))



        #qubit 0 after VCx
        circ.append(transGate, [0])
        circ.append(transGate, [1])
        circ.z(0)
        circ.x(0)
        circ.ry(phi, 0)
        circ.x(0)
        circ.p(15 * np.pi / 16, 0)
        circ.x(0)
        circ.p(3 * np.pi / 16, 0)
        circ.x(0)

        #qubit 1 after VCx
        circ.z(1)
        circ.x(1)
        circ.ry(phi, 1)
        circ.x(1)
        circ.z(1)
        circ.x(1)
        circ.p(13 * np.pi / 16, 1)
        circ.x(1)
        circ.p(9 * np.pi / 16, 1)

        circ.barrier(0)
        circ.barrier(1)

        circ.append(VCxGate.decompose(), [0, 1])


        #qubit 0
        circ.p(15 * np.pi / 16, 0)
        circ.x(0)
        circ.p(3 * np.pi / 16, 0)
        circ.x(0)
        circ.z(0)
        circ.x(0)
        circ.ry(phi, 0)
        circ.x(0)

        #qubit 1
        circ.p(15 * np.pi / 16, 1)
        circ.x(1)
        circ.p(3 * np.pi / 16, 1)
        circ.x(1)
        circ.z(1)
        circ.x(1)
        circ.ry(phi, 1)
        circ.x(1)

        circ.append(transGate.conjugate(), [0])
        circ.append(transGate.conjugate(), [1])

        testWeyly = TwoQubitWeylDecomposition(Operator(circ).to_matrix())
        optiF = testWeyly.circuit()
        optiF.name = "+"

        return optiF

    pi = np.pi
    sqrt = np.sqrt
    log = np.log

    common_log_term = log(-1 / 6 * (-1) ** (1 / 4) * (sqrt(24 - 6 * sqrt(3)) + 1j * (3 + sqrt(3))))
    common_factor = 3j + 1j * sqrt(3) + sqrt(6 * (4 - sqrt(3)))
    denominator = 15 * pi - 4 * 1j * common_log_term

    VCx = np.array([
        [
            -(((-1) ** (5 / 8) * common_factor * (30 * pi - 8 * 1j * common_log_term)) / (24 * denominator)) +
            (3 * (-1) ** (5 / 8) * (-30 * pi + 8 * 1j * common_log_term)) / (2 * common_factor * denominator),
            0,
            0,
            (3 * (-1) ** (5 / 8) * (-30 * pi + 8 * 1j * common_log_term)) / (2 * common_factor * denominator) -
            ((-1) ** (5 / 8) * common_factor * (-30 * pi + 8 * 1j * common_log_term)) / (24 * denominator)
        ],
        [
            0,
            -(1 / 2) * (-1) ** (1 / 8) - (1 / 2) * (-1) ** (5 / 8),
            -(1 / 2) * (-1) ** (1 / 8) + (1 / 2) * (-1) ** (5 / 8),
            0
        ],
        [
            0,
            -(1 / 2) * (-1) ** (1 / 8) + (1 / 2) * (-1) ** (5 / 8),
            -(1 / 2) * (-1) ** (1 / 8) - (1 / 2) * (-1) ** (5 / 8),
            0
        ],
        [
            (3 * (-1) ** (5 / 8) * (-30 * pi + 8 * 1j * common_log_term)) / (2 * common_factor * denominator) -
            ((-1) ** (5 / 8) * common_factor * (-30 * pi + 8 * 1j * common_log_term)) / (24 * denominator),
            0,
            0,
            -(((-1) ** (5 / 8) * common_factor * (30 * pi - 8 * 1j * common_log_term)) / (24 * denominator)) +
            (3 * (-1) ** (5 / 8) * (-30 * pi + 8 * 1j * common_log_term)) / (2 * common_factor * denominator)
        ]
    ], dtype=complex)

    weylDecomCirc = TwoQubitWeylDecomposition(VCx)
    VCxGate = weylDecomCirc.circuit()

    def CZqutrit(transGate):
        circ = QuantumCircuit(4, name="CZ")

        circ.append(transGate, [0])
        circ.append(transGate, [1])
        circ.append(transGate, [2])
        circ.append(transGate, [3])
        omega = 2*np.pi/3

        circ.cp(omega, 0, 2)
        circ.cp(omega, 1, 3)
        circ.cp(-1 * omega, 0, 3)
        circ.cp(-1 * omega, 1, 2)

        circ.append(transGate.conjugate().transpose(), [0])
        circ.append(transGate.conjugate().transpose(), [1])
        circ.append(transGate.conjugate().transpose(), [2])
        circ.append(transGate.conjugate().transpose(), [3])

        return circ

    def createGraphCirc(G, transGate, multi=False):

        n = G.number_of_nodes()
        qubitList = [[2*i, 2*i+1] for i in range(n)]
        test = []
        if multi and hasattr(G, "is_multigraph") and G.is_multigraph():
            for u, v, k in G.edges(keys=True):
                if u != v:
                    test.append([
                        qubitList[u][0],
                        qubitList[u][1],
                        qubitList[v][0],
                        qubitList[v][1]
                    ])
        else:
            for u, v in G.edges():
                if u != v:
                    test.append([
                        qubitList[u][0],
                        qubitList[u][1],
                        qubitList[v][0],
                        qubitList[v][1]
                    ])
        graphCirc = QuantumCircuit(2*n)
        cartanCirc = cartanCircuit(VCxGate, transGate)
        for pair in qubitList:
            graphCirc.append(cartanCirc, pair)
        for edge in test:
            graphCirc.append(CZqutrit(transGate), edge)
        return graphCirc

    basisTransformationMatrix = np.identity(2)
    uniGate = UnitaryGate(basisTransformationMatrix)

    final_circ = createGraphCirc(graph, uniGate)

    return final_circ