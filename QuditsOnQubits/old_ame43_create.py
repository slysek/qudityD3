def create_ame43():
    Fgate = get_f_gate()

    def C3_adder():
        circ = QuantumCircuit(4)

        circ.cx(1,2)
        circ.ccx(1,2,3)
        circ.ccx(1,3,2)
        circ.cx(0,3)
        circ.ccx(0,3,2)
        circ.ccx(0,2,3)
        circ.cz(1,2)
        circ.cz(0,3)

        circ.name = "C3"
        return circ

    C3 = C3_adder()

    ame43circ = QuantumCircuit(8)
    ame43circ.append(Fgate, [4,5])
    ame43circ.append(Fgate, [6,7])
    ame43circ.append(C3, [7, 6, 3, 2])
    ame43circ.append(C3, [7, 6, 1, 0])
    ame43circ.append(C3, [5, 4, 3, 2])
    ame43circ.append(C3, [5, 4, 1, 0])
    ame43circ.append(C3, [5, 4, 1, 0])

    return ame43circ