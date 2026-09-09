from golden_model import encode_bits, pack_bits


def test_known_vector():
    samples = [100, 101, 99, 99, 104, 100, 100, 101]
    bits = encode_bits(samples, k=2)
    assert len(bits) > 16
    assert pack_bits(bits).hex() == "0064de19e6"


def test_extreme_delta_is_supported():
    bits = encode_bits([-32768, 32767], k=2)
    assert bits[:16] == [int(c) for c in "1000000000000000"]
    assert bits[-1] in (0, 1)


if __name__ == "__main__":
    test_known_vector()
    test_extreme_delta_is_supported()
    print("golden self-check: PASS")
