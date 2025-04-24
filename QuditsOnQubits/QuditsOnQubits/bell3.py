import igraph
from qiskit import qpy
from qiskit.circuit import QuantumCircuit
import itertools
from qiskit.primitives import StatevectorSampler
from qiskit.visualization import plot_histogram
import numpy as np

def bell_operator(g: igraph.Graph, prep_circ: QuantumCircuit) -> float:
    d = 3
    N = g.vcount()

    # --- 1. Wszystkie r=1 ---
    for e in g.es:
        e["r"] = 1

    # --- 2. Wybór v1,v2 ---
    v1 = v2 = None
    for e in g.es:
        if e["r"] % d != 0:
            v1, v2 = e.tuple
            break
    nbrs = set(g.neighbors(v1)) - {v1, v2}

    # --- 3. Współczynniki c1[k] (grupa 1) ---
    r12 = g.es[g.get_eid(v1, v2)]["r"] % d
    c1 = {0: 1.0}
    for k in range(1, d):
        Jk = [j for j in nbrs
              if g.es[g.get_eid(v1, j)]["r"] % d == (k * r12) % d]
        c1[k] = 1.0 / (len(Jk) + 1)

    # --- 4. Współczynniki c2[j] (grupa 2) ---
    inv_r12 = pow(r12, -1, d)
    c2 = {}
    for j in nbrs:
        r1j = g.es[g.get_eid(v1, j)]["r"] % d
        kj = (r1j * inv_r12) % d
        Jk = [m for m in nbrs
              if g.es[g.get_eid(v1, m)]["r"] % d == (kj * r12) % d]
        c2[j] = (1 - c1[kj]) / len(Jk)

    # --- 5. Przydział literowych obserwabli ---
    assignment = {i: [] for i in range(N)}
    assignment[v1] = [ "X" + "Z"*k       for k in range(d)         ]
    assignment[v2] = [ "Z"*r12 + "X"*k  for k in range(d)         ]
    for j in nbrs:
        r1j = g.es[g.get_eid(v1, j)]["r"] % d
        assignment[j] = ["Z"*r1j + "X"]
    rest = set(range(N)) - {v1, v2} - nbrs
    for i in rest:
        assignment[i] = ["ZX"]

    # --- 6. Wygeneruj wszystkie obwody pomiarowe i policz <exp> ---
    combos = list(itertools.product(*assignment.values()))
    qmap   = [[2*i,2*i+1] for i in range(N)]
    ghz_circs = []
    for combo in combos:
        meas = QuantumCircuit(2*N)
        for idx, lbl in enumerate(combo):
            gate = to_gates(lbl); gate.name = lbl
            meas.append(gate, qmap[idx])
        qc = prep_circ.copy()
        qc.append(meas, range(2*N))
        qc.measure_all()
        ghz_circs.append(qc)

    sampler = StatevectorSampler()
    bell_exp = []
    for circ in ghz_circs:
        job    = sampler.run([circ], shots=10000)
        counts = job.result()[0].data.meas.get_counts()
        bell_exp.append(np.real(exp_bell(counts, 10000, 0)))

    # --- 7. Budujemy wzory (pattern, weight) według eq. (84) ---
    terms = []

    # grupa 1: G~(1,k)
    for k in range(d):
        pat = [None]*N
        # v1, v2
        pat[v1] = "X" + "Z"*k
        pat[v2] = "Z"*r12 + "X"*k
        # sąsiedzi v1
        for j in nbrs:
            r1j = g.es[g.get_eid(v1, j)]["r"] % d
            r2j = g.es[g.get_eid(v2, j)]["r"] % d
            pat[j] = "Z"*((r1j + k*r2j)%d) + "X"
        # reszta
        for x in rest:
            pat[x] = "ZX"
        terms.append((tuple(pat), c1[k]))

    # grupa 2: G~(2,j)
    for j in nbrs:
        pat = [None]*N
        # v1, v2, j
        r1j = g.es[g.get_eid(v1, j)]["r"] % d
        r2j = g.es[g.get_eid(v2, j)]["r"] % d
        pat[v1] = "Z"*r1j + "X"
        pat[v2] = "Z"*((r12+r2j)%d)
        pat[j]  = "X"
        # pozostałe sąsiedztwa v1\{j}
        for m in nbrs - {j}:
            r1m = g.es[g.get_eid(v1, m)]["r"] % d
            rjm = g.es[g.get_eid(j, m)]["r"] % d
            pat[m] = "Z"*((r1m + rjm)%d) + "X"
        # reszta
        for x in rest:
            pat[x] = "ZX"
        terms.append((tuple(pat), c2[j]))

    # grupa 3: G~(3,k) dla k∈rest
    for k in rest:
        pat = [None]*N
        pat[v1] = "Z"*r12
        pat[v2] = "Z"*(g.es[g.get_eid(v2, k)]["r"]%d)
        for m in nbrs:
            r1m = g.es[g.get_eid(v1, m)]["r"]%d
            pat[m] = "Z"*r1m
        pat[k] = "X"
        for x in rest - {k}:
            pat[x] = "ZX"
        terms.append((tuple(pat), 1.0))

    # --- 8. Zsumuj ---
    weights = [ sum(w for p,w in terms if p==combo) for combo in combos ]
    I_bell  = sum(w*e for w,e in zip(weights, bell_exp))

    return I_bell


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