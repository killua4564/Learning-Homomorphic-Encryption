"""
步驟 5-0d：relinearization —— 把乘完的 3 個 part 壓回 2 個

目的：
  5-0c 每乘一次，密文就多一個 part（3 → 4 → 5），鑰匙也變成 (1, s, s², s³, …)。
  relinearization 用一把「公開的轉換金鑰」把 c2·s² 這一項換成只含 1 和 s 的形式，
  讓密文永遠維持 2 個 part。這一步也會說明「特殊質數 P」到底在做什麼。

方案：
  【轉換金鑰 rlk】本質上是「用 s 加密 P·s²」：
      a ← Z_{P·Q} 均勻亂數，e ← {−3..3}^N
      rlk = (b, a)，b = −a·s + e + P·s²   (mod P·Q)，Q = modulus(len(PRIMES))
      → 檢查：b + a·s = P·s² + e
  【relinearize】密文 (c0, c1, c2)，目前模數 q：
      1. 在模數 P·q 底下算 c2·b 和 c2·a（c2 先中心化）
      2. 兩個結果都「中心化 → 除以 P 四捨五入 → mod q」，得到 d0、d1
      3. 新密文 = (c0 + d0, c1 + d1)，scale、level 不變
  為什麼對：d0 + d1·s ≈ (c2·b + c2·a·s) / P = (c2·P·s² + c2·e) / P = c2·s² + c2·e / P
           → c2·s² 這一項被換成了只用 1 和 s 的形式，多出來的雜訊是 c2·e / P

  為什麼需要 P：c2 的係數大約跟 q 一樣大（約 2^85），如果不乘 P、不除 P，
              多出來的雜訊就是 c2·e，大到蓋掉整個訊息；除以 P 之後變成 c2·e / P，
              只要 P 比 q 大，這一項就小於 1。

  ※ 跟 SEAL 的差別：SEAL 會先把 c2 拆成「mod 每個質數的餘數」分開處理，每一份都 < 2^60，
    所以 P 只要「至少跟單一質數一樣大」（60 bit）。這裡沒有拆，所以 P 要比整個 q 還大，取約 2^90。
  ※ rlk 是公開的，而它活在模數 P·Q 底下，所以 P 也要算進安全上限（HE 標準表的總 bit 數）。

結果：
  - z1 · z2 → relinearize → 只用 (1, s) 解密：結果正確，誤差跟 relinearize 前幾乎一樣
  - P = 1：解出約 1e13 的亂數
    （注意：這裡沿用了大 P 產生的 rlk，壞掉的原因是倍數對不上；
      若要觀察「沒有特殊質數時雜訊 c2·e 爆掉」，P 改成 1 之後要重新呼叫 relin_keygen）
  - z3 連乘「mul → relinearize → rescale」：parts 一直維持 2 個，能乘 3 次，第 4 次 level 用完

觀察：
  - rlk 本身就是一個密文：b + a·s = P·s² + e，也就是「s² 的加密（倍數 P）」。
    relinearize = 用它在密文上算出 c2·s²（密文 × 公開多項式 c2），再 rescale 掉 P
  - rlk 裡的 e 不影響正確性，但影響安全性：rlk 是公開的，沒有 e 時 b ≡ −a·s (mod P)
    是無雜訊的線性方程，直接解出 s
  - 出處：relinearization 來自 Brakerski–Vaikuntanathan 2011（https://eprint.iacr.org/2011/344），
    乘 P 再除 P 的技巧來自 Gentry–Halevi–Smart 2012（https://eprint.iacr.org/2012/099）

執行：python3 step5_0d_relin.py
"""
import itertools
import random

import numpy as np

# ── 5-0a 的編碼 ──────────────────────────────────────────────────────────────
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


# ── 模數鏈 + 特殊質數 ─────────────────────────────────────────────────────────
PRIMES = [33554393, 1048573, 1048571, 1048559]
P = 1237940039285380274899124191   # 特殊質數，約 2^90（比整個 Q ≈ 2^85 還大）


def modulus(level):
    q = 1
    for p in PRIMES[:level]:
        q *= p
    return q


# ── 多項式運算（5-0c）────────────────────────────────────────────────────────
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


# ── 加解密、相乘、rescale（你在 5-0c 寫的）──────────────────────────────────────
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


def decrypt(s, ct):
    parts, scale, level = ct
    q = modulus(level)

    m, si = parts[0], s
    for part in parts[1:]:
        m = poly_add(m, poly_mul(part, si, q), q)
        si = poly_mul(si, s, q)

    return decode(poly_centered(m, q), scale)


def mul(ct1, ct2):
    (parts1, scale1, level1), (parts2, scale2, level2) = ct1, ct2
    level = min(level1, level2)
    q = modulus(level)

    parts = [[] for _ in range(len(parts1) + len(parts2) - 1)]
    for i in range(len(parts1)):
        for j in range(len(parts2)):
            parts[i + j] = poly_add(parts[i + j], poly_mul(parts1[i], parts2[j], q), q)

    return parts, scale1 * scale2, level


def rescale(ct):
    parts, scale, level = ct
    p = PRIMES[level - 1]
    q = modulus(level)

    new_parts = [
        [(v + p // 2) // p for v in poly_centered(part, q)]
        for part in parts
    ]
    return new_parts, scale / p, level - 1


# ── 這一步要寫的 ─────────────────────────────────────────────────────────────
def relin_keygen(s):
    # 在模數 P·Q 底下產生 rlk = (b, a)
    # Q = modulus(len(PRIMES))
    # s² = poly_mul(s, s, P·Q)
    # b = −a·s + e + P·s²   (mod P·Q)
    p, q = P, modulus(len(PRIMES))
    a = [random.randint(0, p * q - 1) for _ in range(N)]
    e = [random.randint(-3, 3) for _ in range(N)]

    neg_a = [-v for v in a]
    p_s2 = [p * v for v in poly_mul(s, s, q)]
    b = poly_add(poly_mul(neg_a, s, p * q), poly_add(e, p_s2, p * q), p * q)
    return (b, a)


def relinearize(ct, rlk):
    # ct 的 parts 必須剛好 3 個 (c0, c1, c2)
    # q = modulus(level)，PQ = P·q
    # rlk 的 b、a 先 mod PQ（rlk 是在最上層產生的，可以直接降下來）
    # c2 中心化後，分別跟 b、a 相乘（mod PQ）
    # 每個係數：中心化 → (v + P // 2) // P → mod q，得到 d0、d1
    # 回傳 ([c0 + d0, c1 + d1], scale, level)
    parts, scale, level = ct
    assert len(parts) == 3
    p, q = P, modulus(level)
    (c0, c1, c2), (b, a) = parts, rlk

    c2 = poly_centered(c2, q)
    def linear_c2(c):
        return [(v + p // 2) // p for v in poly_centered(poly_mul(c, c2, p * q), p * q)]

    c0, c1 = poly_add(c0, linear_c2(b), q), poly_add(c1, linear_c2(a), q)
    return [c0, c1], scale, level


if __name__ == "__main__":
    s = keygen()
    rlk = relin_keygen(s)
    z1 = [1.5, -2.0, 3.25, 0.5]
    z2 = [2.0, 0.5, -1.0, 4.0]
    z3 = [1.1, -0.9, 1.2, 0.8]

    # 2. z1 · z2：relinearize 前後各解密一次，比較誤差
    print(f"{P=}")
    print(f"expected={[i * j for i, j in zip(z1, z2)]}")
    enc_mul = mul(encrypt(s, z1), encrypt(s, z2))
    print(f"沒做 relinearize，decrypt={decrypt(s, enc_mul)}")
    enc_relin = relinearize(enc_mul, rlk)
    print(f"做了 relinearize，decrypt={decrypt(s, enc_relin)}")

    # 3. 把 P 改成 1 重跑第 2 題（可以在檔案上方把 P 暫時改成 1）
    P = 1
    rlk_p1 = relin_keygen(s)
    print(f"\n{P=}")
    print(f"expected={[i * j for i, j in zip(z1, z2)]}")
    enc_mul = mul(encrypt(s, z1), encrypt(s, z2))
    print(f"沒做 relinearize，decrypt={decrypt(s, enc_mul)}")
    enc_relin = relinearize(enc_mul, rlk_p1)
    print(f"做了 relinearize，decrypt={decrypt(s, enc_relin)}")

    # 4. z3 連乘：mul → relinearize → rescale，印出 parts 個數、level、最大誤差
    print("\n連乘實驗")
    P = 1237940039285380274899124191
    ct3, expected = encrypt(s, z3), z3
    for i in range(5):
        decrypted = decrypt(s, ct3)
        max_error = max(abs(x - y) for x, y in zip(expected, decrypted))
        print(f"z3 相乘第 {i} 次")
        print(f"{expected=}")
        print(f"{decrypted=}")
        print(f"len(parts)={len(ct3[0])} level={ct3[2]} 最大誤差 {max_error:.2e}")
        ct3 = rescale(relinearize(mul(ct3, encrypt(s, z3)), rlk))
        expected = [a * b for a, b in zip(expected, z3)]
