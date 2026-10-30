"""
步驟 5-0c：RLWE 密文相乘與 rescale —— 把 4-2、4-3 搬到多項式上

目的：
  一次對 4 個 slot 同時做乘法，並用 rescale 把倍數拉回 Δ。
  順便看到 RLWE 比 LWE 好的地方：相乘後密文「只多一個多項式」，不會像 4-3 那樣長成 3^k。

方案：
  密文表示為 (parts, scale, level)：
    parts ：多項式的 list [c0, c1, c2, ...]，解密時相位 = c0 + c1·s + c2·s² + ...
    scale ：目前的倍數（剛加密是 Δ）
    level ：模數鏈還剩幾個質數（跟 4-3 一樣）

  相乘：把兩個密文都看成「s 的多項式」來乘
      (c0 + c1·s)(d0 + d1·s) = c0·d0 + (c0·d1 + c1·d0)·s + c1·d1·s²
    → 新密文 parts = [c0·d0, c0·d1 + c1·d0, c1·d1]，scale 相乘
    一般情況：新的第 k 個 part = Σ (c_i · d_j)，對所有 i + j = k 加總
    （4-3 的 LWE 版本要做 outer product，長度 3 → 9；這裡因為多項式相乘可以交換順序，
      s 的同次方可以合併，所以 2 → 3。）

  rescale：跟 4-3 一樣，每個 part 的每個係數都「中心化 → 除以 p 四捨五入 → mod 新的模數」，
          scale 除以 p，level 減 1

結果：
  - z1 · z2：相乘後 3 個 part，不 rescale（scale = Δ²）和 rescale 後（scale ≈ Δ）都解得正確
  - 連乘 + rescale：parts 個數 2 → 3 → 4 → 5，每乘一次多 1 個（LWE outer product 是 3 → 9 → 27）
  - z3（|值| ≤ 1.2）：前 3 次正確，第 4 次 level 用完（base 也被除掉），解出 0
  - z1（含 3.25）：第 3 次就全部壞掉 —— 3.25^4 ≈ 111.6，Δ·111.6 ≈ 2^26.8 > base/2 = 2^24

觀察：
  - 新的 parts 個數 = len(parts1) + len(parts2) − 1（把密文看成 s 的多項式，次數相加）
  - rescale 不會多出乘法空間：不 rescale 時是「數字變大」吃掉 q，rescale 時是「模數變小」；
    rescale 的用途是讓 scale 維持在 Δ（才能跟其他密文相加），並讓之後的運算更小更快
  - 只要有一個 slot 溢位，4 個 slot 會全部壞掉：溢位發生在多項式係數上，而每個係數都混著所有 slot
  - rescale 必須用整數運算 (v + p // 2) // p；round(v / p) 會因為 float 只有 53 bit 精度而差上千

執行：python3 step5_0c_mul.py
"""
import itertools
import random

import numpy as np

# ── 5-0a 的編碼（複製過來，讓這個檔案可以獨立執行）─────────────────────────────
N = 8
DELTA = 2**20
ZETA = np.exp(1j * np.pi / N)

SLOT_ROOTS = np.array([ZETA ** pow(5, j, 2 * N) for j in range(N // 2)])
SLOT_ROOTS_MATRIX = np.array([SLOT_ROOTS ** i for i in range(N)]).transpose()

ALL_ROOTS = np.concatenate([SLOT_ROOTS, np.conj(SLOT_ROOTS)])
ALL_ROOTS_MATRIX = np.array([ALL_ROOTS ** i for i in range(N)]).transpose()


def encode(z, delta):
    target = np.concatenate([z, np.conj(z)])
    m = np.linalg.solve(ALL_ROOTS_MATRIX, target)
    return [int(v) for v in (m.real * delta).round()]


def decode(m, delta):
    m = np.array(m, dtype=np.float64).reshape(N, 1)
    z = np.matmul(SLOT_ROOTS_MATRIX, m)
    return [float(v) for v in (z.real / delta).reshape(N // 2)]


# ── 模數鏈（跟 4-3 一樣：base ≈ 2^25，後面三個質數離 2^20 只差 3、5、17）──────────
PRIMES = [33554393, 1048573, 1048571, 1048559]


def modulus(level):
    q = 1
    for p in PRIMES[:level]:
        q *= p
    return q


# ── 5-0b 的多項式運算（加上 q 參數，因為模數會隨 level 改變）─────────────────────
def poly_add(a, b, q):
    return [(x + y) % q for x, y in itertools.zip_longest(a, b, fillvalue=0)]


def poly_mul(a, b, q):
    c = [0] * N
    for i in range(N):
        for j in range(N):
            if i + j < N:
                c[i + j] += a[i] * b[j]
            else:
                c[i + j - N] -= a[i] * b[j]
    return [v % q for v in c]


def poly_centered(a, q):
    return [(v + q // 2) % q - q // 2 for v in a]


# ── 5-0b 的加密（改成回傳 (parts, scale, level)）────────────────────────────────
def keygen():
    return [random.randint(-1, 1) for _ in range(N)]


def encrypt(s, z):
    level = len(PRIMES)
    q = modulus(level)
    m = encode(z, DELTA)
    a = [random.randint(0, q - 1) for _ in range(N)]
    e = [random.randint(-3, 3) for _ in range(N)]

    neg_a = [-v for v in a]
    c = poly_add(poly_mul(neg_a, s, q), poly_add(m, e, q), q)
    return [c, a], DELTA, level


# ── 這一步要寫的 ─────────────────────────────────────────────────────────────
def decrypt(s, ct):
    # parts, scale, level = ct
    # 相位 = parts[0] + parts[1]·s + parts[2]·s² + ...（都在模數 modulus(level) 底下）
    # 提示：用一個變數記住「s 的 i 次方」，每處理一個 part 就再乘一次 s
    # 中心化後 decode(相位, scale)
    parts, scale, level = ct
    q = modulus(level)

    m, si = parts[0], s
    for part in parts[1:]:
        m = poly_add(m, poly_mul(part, si, q), q)
        si = poly_mul(si, s, q)

    return decode(poly_centered(m, q), scale)


def mul(ct1, ct2):
    # 兩邊 level 不同時，取較低的 level，把較高那邊的係數 mod 較小的模數
    # 新的第 k 個 part = Σ poly_mul(c[i], d[j])，對所有 i + j = k 加總
    # scale 相乘
    (parts1, scale1, level1), (parts2, scale2, level2) = ct1, ct2
    level = min(level1, level2)
    q = modulus(level)

    parts = [[] for _ in range(len(parts1) + len(parts2) - 1)]
    for i in range(len(parts1)):
        for j in range(len(parts2)):
            parts[i + j] = poly_add(parts[i + j], poly_mul(parts1[i], parts2[j], q), q)

    return parts, scale1 * scale2, level


def rescale(ct):
    # p = PRIMES[level − 1]，每個 part 的每個係數：中心化 → 除以 p 四捨五入 → mod modulus(level − 1)
    # scale / p，level − 1
    # 不能寫 round(v / p)：v 大約 2^85，v / p 會變成 float，只有 53 bit 精度，
    # 四捨五入會差上千（4-3 用 Sage 的有理數所以沒事）。要用整數運算：(v + p // 2) // p
    parts, scale, level = ct
    p = PRIMES[level - 1]
    q = modulus(level)

    new_parts = [
        [(v + p // 2) // p for v in poly_centered(part, q)]
        for part in parts
    ]
    return new_parts, scale / p, level - 1


if __name__ == "__main__":
    s = keygen()
    z1 = [1.5, -2.0, 3.25, 0.5]
    z2 = [2.0, 0.5, -1.0, 4.0]

    # 2. z1 · z2
    # 不 rescale 解密、rescale 後解密，印出結果、parts 個數、level
    expected = [a * b for a, b in zip(z1, z2)]
    enc_mul = mul(encrypt(s, z1), encrypt(s, z2))
    enc_rescale = rescale(enc_mul)
    print(f"{expected=}")
    print(f"不 rescale 解密 {decrypt(s, enc_mul)}，len(parts)={len(enc_mul[0])}，level={enc_mul[2]}")
    print(f"有 rescale 解密 {decrypt(s, enc_rescale)}，len(parts)={len(enc_rescale[0])}，level={enc_rescale[2]}")

    # 3. 連乘
    z3 = [1.1, -0.9, 1.2, 0.8]
    # 對 z1 和 z3 各做一次：ct = encrypt(s, z)，連續 ct = rescale(mul(ct, encrypt(s, z)))，
    # 印出 parts 個數、level、和正確值的最大誤差；level 用完時會怎樣？
    ct1, expected = encrypt(s, z1), z1
    for i in range(1, 5):
        ct1 = rescale(mul(ct1, encrypt(s, z1)))
        expected = [a * b for a, b in zip(expected, z1)]
        print(f"z1 相乘第 {i} 次")
        print(f"{expected=}")
        print(f"decrypt={decrypt(s, ct1)}")

    ct3, expected = encrypt(s, z3), z3
    for i in range(1, 5):
        ct3 = rescale(mul(ct3, encrypt(s, z3)))
        expected = [a * b for a, b in zip(expected, z3)]
        print(f"z3 相乘第 {i} 次")
        print(f"{expected=}")
        print(f"decrypt={decrypt(s, ct3)}")
