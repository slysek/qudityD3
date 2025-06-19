from qiskit import qpy

def create_ame43():

    with open('ame43fig12.qpy', 'rb') as fd:
        ame43circ = qpy.load(fd)[0]

    return ame43circ