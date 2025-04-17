from qiskit_aer import AerSimulator
from qiskit_ibm_runtime.fake_provider import FakeBrisbane
from qiskit.transpiler.preset_passmanagers import generate_preset_pass_manager
from qiskit_ibm_runtime import SamplerV2

def brisbane(circ):
    backend_brisbane = FakeBrisbane()

    target_brisbane = backend_brisbane.target
    pass_manager_brisbane = generate_preset_pass_manager(
        optimization_level=3, target=target_brisbane
    )

    tempCirc = circ.copy()
    tempCirc.measure_all()
    transpiled_brisbane = pass_manager_brisbane.run(tempCirc)

    transpiled_brisbane.name = "brisbane"
    # aer = AerSimulator.from_backend(backend_brisbane)
    # sampler = SamplerV2(aer)
    # job = sampler.run([transpiled_brisbane])
    # pub_result = job.result()[0]
    # counts = pub_result.data.meas.get_counts()

    return transpiled_brisbane