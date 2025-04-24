import numpy as np
import igraph
from igraph import Graph
from qiskit import qpy
from qiskit.circuit import QuantumCircuit
from qiskit.quantum_info import Statevector, Operator
from qiskit.circuit.library import UnitaryGate
from qiskit.primitives import StatevectorSampler
import itertools

def get_f_gate():
    with open('Fgate.qpy', 'rb') as fd:
        return qpy.load(fd)[0]

def create_ghz(raw=False):
    with open('GHZ.qpy', 'rb') as fd:
        ghz_circ = qpy.load(fd)[0]
    fgate = get_f_gate()
    if raw:
        return ghz_circ, Statevector(ghz_circ)
    ghz_prep = ghz_circ.copy()
    ghz_prep.remove_final_measurements()
    cartanH = Operator(fgate).conjugate().transpose()
    change = cartanH.tensor(cartanH).tensor(Operator(np.identity(4))).data
    sv = Statevector(ghz_prep)
    vec = sv.data.reshape((-1,1))
    new_vec = change @ vec
    return ghz_circ, Statevector(new_vec.flatten())

def create_graph_ghz(graph: Graph) -> QuantumCircuit:
    with open('Fgate.qpy', 'rb') as fd:
        cartanCirc = qpy.load(fd)[0]
    with open('CZ.qpy', 'rb') as fd:
        CZ = qpy.load(fd)[0]
    def createGraphCirc(G, cartan):
        n = G.vcount()
        qc = QuantumCircuit(2*n)
        for i in range(n):
            qc.append(cartan, [2*i, 2*i+1])
        for u,v in G.get_edgelist():
            qc.append(CZ, [2*u,2*u+1,2*v,2*v+1])
        return qc
    return createGraphCirc(graph, cartanCirc)

def to_gates(s: str) -> QuantumCircuit:
    with open('Zgate.qpy', 'rb') as fd:
        Zgate = qpy.load(fd)[0]
    with open('Xgate.qpy', 'rb') as fd:
        Xgate = qpy.load(fd)[0]
    qc = QuantumCircuit(2)
    for ch in s:
        if ch == 'X':
            qc.append(Xgate, [0,1])
        elif ch == 'Z':
            qc.append(Zgate, [0,1])
    return qc

def exp_bell_statevector(statevec: np.ndarray) -> float:
    d = 3
    omega = np.exp(2j*np.pi/d)
    expv = 0+0j
    n_qubits = int(np.log2(statevec.size))
    for idx, amp in enumerate(statevec):
        p = abs(amp)**2
        if p==0: continue
        bits = format(idx, f'0{n_qubits}b')
        s = 0
        for i in range(0, n_qubits, 2):
            pair = bits[i:i+2]
            if pair=='00': v=0
            elif pair=='01': v=1
            elif pair=='10': v=2
            else: v=0
            s += v
        expv += p * omega**s
    return expv.real

def bell_operator(g: Graph, circuit: QuantumCircuit):
    d = 3
    N = g.vcount()
    # przypisz r=1
    for e in g.es: e['r']=1
    # wybierz v1,v2
    v1=v2=None
    for e in g.es:
        if e['r']%d:
            v1,v2=e.tuple
            break
    nbrs = set(g.neighbors(v1)) - {v1,v2}
    # oblicz c_{i,k}
    r12 = g.es[g.get_eid(v1,v2)]['r']%d
    coeffs={}
    for k in range(1,d):
        Jk=[j for j in nbrs if g.es[g.get_eid(v1,j)]['r']%d == (k*r12)%d]
        m=len(Jk)
        c=1/(m+1)
        coeffs[(v1,k)] = c
        for j in Jk: coeffs[(j,k)] = c
    print("Współczynniki c_{i,k}:")
    for (v,k),c in sorted(coeffs.items()):
        print(f" c_{{{v+1},{k}}} = {c:.3f}")
    # przydział obserwabli
    assign={i:[] for i in range(N)}
    assign[v1]=['X'+'Z'*k for k in range(d)]
    assign[v2]=['Z'*r12+'X'*k for k in range(d)]
    for i in nbrs:
        r = g.es[g.get_eid(v1,i)]['r']%d
        assign[i]=['Z'*r+'X']
    for i in set(range(N)) - {v1,v2} - nbrs:
        assign[i]=['ZX']
    # obwiednie rotacji (bez pomiarów)
    meas_circs=[]
    qlist=[[2*i,2*i+1] for i in range(N)]
    for combo in itertools.product(*assign.values()):
        mc=QuantumCircuit(2*N)
        for idx,gs in enumerate(combo):
            gate=to_gates(gs); gate.name=gs
            mc.append(gate, qlist[idx])
        meas_circs.append(mc)
    # sampler + analityczne <Bell>
    sampler=StatevectorSampler()
    vals=[]
    for mc in meas_circs:
        full = circuit.copy()
        full.remove_final_measurements()      # <--- teraz w miejscu
        full.append(mc, range(2*N))
        job = sampler.run([full])
        sv = job.result().get_statevector(0)
        vals.append(exp_bell_statevector(sv))
    return vals, meas_circs
