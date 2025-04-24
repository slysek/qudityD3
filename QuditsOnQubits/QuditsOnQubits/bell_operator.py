import igraph
from qiskit import qpy
from qiskit.circuit import QuantumCircuit
import itertools
from qiskit.primitives import StatevectorSampler
from qiskit.visualization import plot_histogram
import numpy as np

def bell_operator(g: igraph.Graph, circuit: QuantumCircuit):
    """
    Algorytm:
      1. Znajduje dwa wierzchołki (v1,v2) połączone krawędzią r_{v1,v2} != 0 (mod 3).
      2. Określa sąsiedztwo v1 (bez v2).
      3. Przydziela:
         - v1: obserwable typu XZ^k (k=0..2)
         - v2: obserwable typu Z^(r_{v1,v2}) X^k (k=0..2)
         - inni sąsiedzi v1 (poza v2): np. ["Z^r_{v1,i} X^1"]
         - pozostałe wierzchołki: np. ["Z^1", "X^1"]

    """

    d = 3
    N = g.vcount()



    for e in g.es:
        e["r"] = 1

    # KROK A: wybór (v1, v2) - pierwsza krawędź z r != 0.
    v1, v2 = None, None
    for e in g.es:
        if (e["r"] % d) != 0:
            v1, v2 = e.source, e.target
            break

    # KROK B: sąsiedztwo v1 (oprócz v2)
    neighbors_v1 = set(g.neighbors(v1, mode="all"))
    if v2 in neighbors_v1:
        neighbors_v1.remove(v2)
    if v1 in neighbors_v1:
        neighbors_v1.remove(v1)

    # KROK C: słownik wynikowy: {wierzchołek: [lista_obserwabli]}
    assignment = { i: [] for i in range(N) }

    # KROK D: v1 -> XZ^k, k=0..2
    obs_v1 = []
    for k in range(d):

        str_v1 = "X"
        for j in range(k):
            str_v1 += "Z"
        obs_v1.append(str_v1)

    assignment[v1] = obs_v1

    # KROK E: v2 -> Z^(r_{v1,v2}) X^k, k=0..2
    r_v1v2 = g.es[g.get_eid(v1, v2)]["r"] % 3
    obs_v2 = []
    for k in range(d):
        str_v2 = ""

        for j in range(r_v1v2):
            str_v2 += "Z"

        for j2 in range(k):
            str_v2 += "X"

        obs_v2.append(str_v2)

    assignment[v2] = obs_v2

    # KROK F: inni sąsiedzi v1: "Z^r_{v1,i} X" (jeden operator lub więcej)

    for i in neighbors_v1:
        r_v1i = g.es[g.get_eid(v1, i)]["r"] % 3
        str_neighbors = ""

        for j in range(r_v1i):
            str_neighbors += "Z"
        str_neighbors += "X"

        assignment[i] = [str_neighbors]

    # KROK G: reszta wierzchołków
    all_verts = set(range(N))
    remainder = all_verts - {v1, v2} - neighbors_v1

    circ_reminder_list = []
    for i in remainder:
        assignment[i] = ["ZX"]  # dwie obserwable bazowe




    final_circ_list = []

    qutrit_list = [[2 * i, 2 * i + 1] for i in range(N)]

    str_list_gates = []

    for combination in itertools.product(*list(assignment.values())):
        circ = QuantumCircuit(2*N)
        for j in range(len(combination)):
            str_list_gates.append(combination[j])
            tempcirc = to_gates(combination[j])
            tempcirc.name = combination[j]
            circ.append(tempcirc, qutrit_list[j])
        final_circ_list.append(circ)

    #print(assignment)

    ghz_bell_circs = []

    for i in final_circ_list:
        tempcirc2 = circuit.copy()
        tempcirc2.append(i, [k for k in range(2*N)])
        tempcirc2.measure_all()
        ghz_bell_circs.append(tempcirc2)

    shots = 10000
    sampler = StatevectorSampler()

    bell_exp_list = []
    exps_operators = []
    for i, j in zip(ghz_bell_circs, range(len(ghz_bell_circs))):
        job = sampler.run([i], shots = shots)
        data_pub = job.result()[0].data
        counts = data_pub.meas.get_counts()
        exps_operators.append(str_list_gates[j])
        bell_exp = exp_bell(counts, shots, j)
        bell_exp_list.append(np.real(bell_exp))

    print(assignment)
    print(dict(zip(exps_operators, bell_exp_list)))
    return bell_exp_list, ghz_bell_circs

def to_gates(str):

    with open('Zgate.qpy', 'rb') as fd:
        Zgate = qpy.load(fd)[0]

    with open('Xgate.qpy', 'rb') as fd:
        Xgate = qpy.load(fd)[0]

    circ = QuantumCircuit(2)
    for i in str:
        if i == "X":
            circ.append(Xgate, [0, 1])
        if i == "Z":
            circ.append(Zgate, [0, 1])

    return circ

def exp_bell(counts, shots, k):
    oper = 0
    keys = list(counts.keys())[k]
    for i in list(counts.keys()):
        p = counts[i]/shots
        sum = 0
        for j in range(0, len(i), 2):
            if f'{keys[j]}{keys[j+1]}' == "00":
                sum += 0
            elif f'{keys[j]}{keys[j+1]}' == "01":
                sum += 1
            elif f'{keys[j]}{keys[j+1]}' == "10":
                sum += 2
        omega = np.exp(2*np.pi*1.j/3)
        f = pow(omega, sum)
        oper += p*f

    return oper