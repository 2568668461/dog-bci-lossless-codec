"""Bit-exact reference for rice_demo_encoder.v."""

from __future__ import annotations

from typing import Iterable, List


def zigzag17(value: int) -> int:
    if not -(1 << 16) <= value <= (1 << 16) - 1:
        raise ValueError(f"delta outside signed 17-bit range: {value}")
    return ((value << 1) ^ (value >> 16)) & ((1 << 17) - 1)


def encode_bits(samples: Iterable[int], k: int = 2) -> List[int]:
    values = [int(v) for v in samples]
    if not values:
        raise ValueError("at least one sample is required")
    if any(v < -32768 or v > 32767 for v in values):
        raise ValueError("samples must be signed 16-bit values")
    if not 0 <= k <= 16:
        raise ValueError("k must be in [0, 16]")

    bits = [(values[0] >> bit) & 1 for bit in range(15, -1, -1)]
    for previous, current in zip(values, values[1:]):
        code = zigzag17(current - previous)
        quotient, remainder = code >> k, code & ((1 << k) - 1 if k else 0)
        bits.extend([0] * quotient)
        bits.append(1)
        bits.extend((remainder >> bit) & 1 for bit in range(k - 1, -1, -1))
    return bits


def pack_bits(bits: Iterable[int]) -> bytes:
    result = bytearray()
    acc = 0
    count = 0
    for bit in bits:
        acc = (acc << 1) | int(bit)
        count += 1
        if count == 8:
            result.append(acc)
            acc = count = 0
    if count:
        result.append(acc << (8 - count))
    return bytes(result)


if __name__ == "__main__":
    vector = [100, 101, 99, 99, 104, 100, 100, 101]
    bits = encode_bits(vector)
    print("samples:", vector)
    print("bits:", "".join(map(str, bits)))
    print("packed:", pack_bits(bits).hex(" "))
