import igraph
from qiskit import qpy
from qiskit.circuit import QuantumCircuit
from qiskit.primitives import StatevectorSampler
import numpy as np


def to_gates(sym: str) -> QuantumCircuit:
    """
    Ładuje bramkę X lub Z^k z pliku QPY do jednego qutrytu (zakodowanego na 2 qubitach).
    """
    with open('Zgate.qpy', 'rb') as fd:
        Zgate = qpy.load(fd)[0]
    with open('Xgate.qpy', 'rb') as fd:
        Xgate = qpy.load(fd)[0]

    circ = QuantumCircuit(2)
    for ch in sym:
        if ch == 'X':
            circ.append(Xgate, [0, 1])
        elif ch == 'Z':
            circ.append(Zgate, [0, 1])
    return circ


def exp_bell(counts: dict, shots: int) -> complex:
    """
    Oblicza wartość oczekiwaną exp_bell z histogramu counts dla qutrytów.
    Ignoruje nielegalne pary '11'.
    """
    ω = np.exp(2j * np.pi / 3)
    mapping = {'00': 0, '01': 1, '10': 2}
    total = sum(counts.values())
    oper = 0+0j
    for bitstr, cnt in counts.items():
        s = 0
        valid = True
        for i in range(0, len(bitstr), 2):
            pair = bitstr[i:i+2]
            if pair not in mapping:
                valid = False
                break
            s += mapping[pair]
        if not valid:
            continue
        p = cnt / total
        oper += p * (ω ** s)
    return oper


def bell_operator(g: igraph.Graph, circuit: QuantumCircuit):
    d = 3
    N = g.vcount()

    # --- 1) Ustaw wagę r=1 na każdej krawędzi (mod d) ---
    if 'r' not in g.es.attribute_names():
        g.es['r'] = [1] * g.ecount()
    g.es['r'] = [int(r) % d for r in g.es['r']]

    # --- 2) Wybór (v1, v2) pierwszej krawędzi z r != 0 ---
    v1 = v2 = None
    for e in g.es:
        if e['r'] != 0:
            v1, v2 = e.tuple
            break

    # --- 3) Sąsiedztwo v1 (bez v2) ---
    neighbors_v1 = set(g.neighbors(v1, mode='all'))
    neighbors_v1.discard(v1)
    neighbors_v1.discard(v2)

    # --- 4) Tworzenie assignment jak dotychczas ---
    assignment = {i: [] for i in range(N)}
    # D: v1 -> X Z^k
    obs_v1 = ['X' + 'Z'*k for k in range(d)]
    assignment[v1] = obs_v1
    # E: v2 -> Z^r X^k
    r12 = g.es[g.get_eid(v1, v2)]['r']
    obs_v2 = ['Z'*r12 + 'X'*k for k in range(d)]
    assignment[v2] = obs_v2
    # F: inni sąsiedzi v1
    for j in neighbors_v1:
        r1j = g.es[g.get_eid(v1, j)]['r']
        assignment[j] = ['Z'*r1j + 'X']
    # G: reszta
    for i in set(range(N)) - {v1, v2} - neighbors_v1:
        assignment[i] = ['ZX']

    # --- 5) Pełna tomografia: budowa i uruchomienie obwodów ---
    qutrit_list = [[2*i, 2*i+1] for i in range(N)]
    circuits = []
    str_list_gates = []
    import itertools
    for combo in itertools.product(*assignment.values()):
        circ = circuit.copy()
        for idx, gate_str in enumerate(combo):
            str_list_gates.append(gate_str)
            temp = to_gates(gate_str)
            circ.append(temp, qutrit_list[idx])
        circ.measure_all()
        circuits.append(circ)

    sampler = StatevectorSampler()
    job = sampler.run(circuits)
    results = job.result()

    bell_exp_list = []
    exps_operators = []
    for i, gate_str in enumerate(str_list_gates):
        counts = results[i].data.meas.get_counts()
        exps_operators.append(gate_str)
        bell_exp_list.append(np.real(exp_bell(counts, shots=10000)))

    operators = ['--'.join(exps_operators[i:i+N]) for i in range(0, len(exps_operators), N)]
    print('Wartości oczekiwane:')
    for op, val in zip(operators, bell_exp_list):
        print(f" {op} -> {val:.6f}")

    # --- 6) Obliczanie współczynników c_{i,k} ---
    coeffs = {}
    r_ref = g.es[g.get_eid(v1, v2)]['r']
    for k in range(1, d):
        Jk = [j for j in neighbors_v1
              if g.es[g.get_eid(v1, j)]['r'] == (k * r_ref) % d]
        c_val = 1.0 / (len(Jk) + 1)
        coeffs[(v1, k)] = c_val
        for j in Jk:
            coeffs[(j, k)] = c_val

    print('Współczynniki c_{i,k}:')
    for (v, k), c in coeffs.items():
        print(f" c_{{{v},{k}}} = {c:.3f}")

    # --- 7) Obliczenie finalnego I_Bell ---
    exp_map = dict(zip(operators, bell_exp_list))
    I_bell = 0.0
    for (vert, k), c in coeffs.items():
        combo = []
        for j in range(N):
            if j == vert:
                combo.append(assignment[j][k])
            else:
                combo.append(assignment[j][0])
        key = '--'.join(combo)
        E = exp_map.get(key, 0.0)
        print(f"Dla G_{{{vert},{k}}}: key={key}, c={c:.3f}, E={E:.6f}")
        I_bell += c * E

    print(f"Finalne I_Bell = {I_bell:.6f}")
    return bell_exp_list, circuits, I_bell
