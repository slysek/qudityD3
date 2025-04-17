from qiskit_aer import AerSimulator
from qiskit_ibm_runtime.fake_provider import FakeTorino
from qiskit.transpiler.preset_passmanagers import generate_preset_pass_manager
from qiskit_ibm_runtime import SamplerV2

def torino(circ):

    backend_torino = FakeTorino()

    target_torino = backend_torino.target
    pass_manager_torino = generate_preset_pass_manager(
        optimization_level=3, target=target_torino
    )
    tempCirc = circ.copy()
    tempCirc.measure_all()
    transpiled_torino = pass_manager_torino.run(tempCirc)

    transpiled_torino.name = "torino"

    # aer = AerSimulator.from_backend(backend_torino)
    # sampler = SamplerV2(aer)
    # job = sampler.run([transpiled_torino])
    # pub_result = job.result()[0]
    # counts = pub_result.data.meas.get_counts()
    #
    # self.statTorino = counts

    return transpiled_torino