"""Tests for the repetition-code QEC module."""

import numpy as np
import pytest

from qec_vqe.qec.repetition_code import RepetitionCode, memory_experiment


@pytest.mark.parametrize("bad", [(2, "X"), (4, "X"), (3, "H")])
def test_invalid_code_specs(bad):
    d, basis = bad
    with pytest.raises(ValueError):
        RepetitionCode(d, basis=basis)


def test_layout_counts():
    c = RepetitionCode(5, basis="X")
    assert c.num_qubits == 9  # 2d-1
    assert c.num_data == 5
    assert c.num_ancilla == 4
    assert c.data_qubits() == [0, 2, 4, 6, 8]
    assert c.ancilla_qubits() == [1, 3, 5, 7]


def test_decode_single_error():
    c = RepetitionCode(5, basis="X")
    # single error at data index k -> syndrome bits (k-1, k); data bit k flips.
    # error at data index 1: syndrome [1, 1, 0, 0], flipped data bit 1
    assert c.decode([1, 1, 0, 0], [0, 1, 0, 0, 0]) == 0
    # error at data index 2: syndrome [0, 1, 1, 0], flipped data bit 2
    assert c.decode([0, 1, 1, 0], [0, 0, 1, 0, 0]) == 0
    # no error: all-zero syndrome, majority vote of unflipped data
    assert c.decode([0, 0, 0, 0], [0, 0, 0, 0, 0]) == 0
    assert c.decode([0, 0, 0, 0], [1, 1, 1, 1, 1]) == 1


def test_noiseless_memory_no_errors():
    """Without noise the code must return zero logical error."""
    for basis in ("X", "Z"):
        c = RepetitionCode(3, basis=basis)
        r = memory_experiment(c, p_x=0.0, p_y=0.0, p_z=0.0, shots=1000, seed=1)
        assert r["logical_error_rate"] == 0.0


def test_xcode_suppresses_z_errors_with_distance():
    """pL must fall monotonically with distance under pure Z (phase) noise."""
    rates = []
    for d in (3, 5, 7):
        c = RepetitionCode(d, basis="X")
        r = memory_experiment(c, p_x=0.0, p_y=0.0, p_z=0.05, rounds=1,
                              shots=100000, seed=3)
        rates.append(r["logical_error_rate"])
    # rates strictly decrease
    assert rates[0] > rates[1] > rates[2]
    # d=3 theoretical value is 3 p_z^2 = 0.0075 (within binomial error)
    se = np.sqrt(rates[0] * (1 - rates[0]) / 100000)
    assert abs(rates[0] - 3 * 0.05 ** 2) < 6 * se


def test_zcode_suppresses_x_errors_with_distance():
    r3 = memory_experiment(RepetitionCode(3, basis="Z"), p_x=0.05, p_y=0.0,
                           p_z=0.0, rounds=1, shots=60000, seed=3)
    r5 = memory_experiment(RepetitionCode(5, basis="Z"), p_x=0.05, p_y=0.0,
                           p_z=0.0, rounds=1, shots=60000, seed=3)
    assert r5["logical_error_rate"] < r3["logical_error_rate"]


def test_error_grows_with_rounds():
    """More syndrome rounds at fixed physical rate = more exposure time."""
    c = RepetitionCode(3, basis="X")
    one = memory_experiment(c, p_x=0.001, p_y=0.001, p_z=0.02, rounds=1,
                            shots=30000, seed=5)
    four = memory_experiment(c, p_x=0.001, p_y=0.001, p_z=0.02, rounds=4,
                             shots=30000, seed=5)
    assert four["logical_error_rate"] > one["logical_error_rate"]


def test_gate_overhead_linear_in_distance():
    """KPI: physical gates per syndrome cycle scale linearly (sub-quadratic)."""
    prev = 0
    for d in (3, 5, 7, 9):
        c = RepetitionCode(d, basis="X")
        g = c.physical_gates_per_cycle()
        assert g > prev  # monotonic
        # linear in d: gates ~ (a*d + b); check ratio bounds
        assert g < 20 * d  # very loose linear bound (empirical ~6d+10)
        prev = g
    # explicit sub-quadratic check: gates <= c * qubits^2 is trivially true,
    # so instead verify gates/qubit stays roughly constant (linear overhead)
    ratios = [
        RepetitionCode(d, basis="X").physical_gates_per_cycle() /
        RepetitionCode(d, basis="X").num_qubits
        for d in (3, 5, 7, 9)
    ]
    assert max(ratios) / min(ratios) < 2.0  # ~flat => linear scaling


def test_cx_count_scales_linearly():
    cx = [RepetitionCode(d, basis="X").gates_per_cycle()["cx"] for d in (3, 5, 7, 9)]
    assert cx == [4, 8, 12, 16]
