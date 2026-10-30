"""
步驟 5-0b：RLWE 加解密 —— 一個密文裝 4 個數字

目的：
  把 5-0a 編碼好的多項式加密。方案跟 4-1 一模一樣，只是「數字」換成「mod X^N + 1 的多項式」：

      4-1（LWE，一個數字）                       5-0b（RLWE，一個多項式 = 4 個 slot）
      s ← {−1, 0, 1}^n                         s(X)：N 個係數各取 {−1, 0, 1}
      a ← Z_q^n                                a(X)：N 個係數各取 Z_q 的均勻亂數
      e ← {−3, …, 3}                           e(X)：N 個係數各取 {−3, …, 3}
      c = (−<a,s> + round(Δx) + e, a)          c = (−a·s + m + e, a)，m = encode(z, Δ)
      相位 = c0 + <a, s>                       相位 = c0 + c1·s  （多項式乘法 mod X^N + 1，係數 mod q）

  密文永遠是 2 個多項式 (c0, c1)，不管裡面裝幾個數字。
  解密時的「相位」c0 + c1·s，就是 4-1 的 c0 + <a, s>：內積換成多項式乘法而已。

  為什麼只有乘法要處理 X^N + 1：slot 點都是 X^N + 1 的根（ζ^N = −1），
  在這些點上 X^N 和 −1 的值一樣，所以把 X^N 換成 −1 不會改變 slot 的值。

結果：
  - 加密 z1 再解密，以及兩個密文相加，誤差約 1e-6
  - 每個 slot 的平均誤差 3.15e-6，乘回 Δ 約 3.3；4-1 單一數字時約 1.7（E|e| = 12/7）
  - sample extraction：只重新排列密文係數（a_vec = [c1[0], −c1[N−1], …, −c1[1]]，即負循環矩陣第一列），
    就得到 slot 總和 × 2/N 的 LWE 密文，解出的總和誤差約 1e-5，完全不需要 galois key

觀察：
  - 加密用的乘法必須是多項式乘法（mod X^N + 1）。若誤用逐係數相乘，加解密仍會「看起來正確」，
    但等於 8 個各自獨立的 1 維 LWE，知道一組明文就能逐格試出 s，而且之後的密文相乘會壞掉
  - 每個 slot 的雜訊是 N 個係數雜訊的加權總和，所以比單一數字的 LWE 大（詳細推導列在待補清單）
  - sample extraction 的結果是單一數字的 LWE 密文，已離開「一個密文裝很多數字」的形式，
    適合當最後輸出；需要繼續逐格運算時，SEAL 改用旋轉（5-0e）

執行：python3 step5_0b_rlwe.py
"""
import random

import numpy as np

# ── 以下是 5-0a 寫好的編碼（直接複製過來，讓這個檔案可以獨立執行）──────────────
N = 8
DELTA = 2**20
ZETA = np.exp(1j * np.pi / N)

SLOT_ROOTS = np.array([ZETA ** pow(5, j, 2 * N) for j in range(N // 2)])
SLOT_ROOTS_MATRIX = np.array([SLOT_ROOTS ** i for i in range(N)]).transpose()

ALL_ROOTS = np.concatenate([SLOT_ROOTS, np.conj(SLOT_ROOTS)])
ALL_ROOTS_MATRIX = np.array([ALL_ROOTS ** i for i in range(N)]).transpose()


def add(a, b):
    return [i + j for i, j in zip(a, b)]


def mul(a, b):
    return [i * j for i, j in zip(a, b)]


def encode(z, delta):
    # 4 個 slot 值 → N 個整數係數（回傳 Python int 的 list，避免之後跟 2^60 的數運算時溢位）
    target = np.concatenate([z, np.conj(z)])
    m = np.linalg.solve(ALL_ROOTS_MATRIX, target)
    return [int(v) for v in (m.real * delta).round()]


def decode(m, delta):
    # N 個係數 → 代入 4 個 slot 點，除以 delta
    m = np.array(m, dtype=np.float64).reshape(N, 1)
    z = np.matmul(SLOT_ROOTS_MATRIX, m)
    return [float(v) for v in (z.real / delta).reshape(N // 2)]
# ─────────────────────────────────────────────────────────────────────────────

Q = 2**60 + 33   # 約 2^60 的奇數（這一步只有一個模數，還沒有模數鏈）


def poly_add(a, b, q=Q):
    return [(x + y) % q for x, y in zip(a, b)]


def poly_mul(a, b, q=Q):
    # 多項式相乘 mod (X^N + 1)，係數 mod q
    c = [0] * N
    for i in range(N):
        for j in range(N):
            if i + j < N:
                c[i + j] += a[i] * b[j]
            else:
                c[i + j - N] -= a[i] * b[j]
    return [v % q for v in c]


def poly_centered(a, q=Q):
    return [(v + q // 2) % q - q // 2 for v in a]


def keygen():
    # 回傳 s，長度 N 的 list，每個係數取自 {−1, 0, 1}
    return [random.randint(-1, 1) for _ in range(N)]


def encrypt(s, z):
    # m = encode(z, DELTA)
    # a = N 個 Z_q 均勻亂數，e = N 個 {−3..3}
    # c0 = −a·s + m + e（逐係數 mod Q），c1 = a
    # 回傳 (c0, c1)
    m = encode(z, DELTA)
    a = [random.randint(0, Q - 1) for _ in range(N)]
    e = [random.randint(-3, 3) for _ in range(N)]

    neg_a = [-v for v in a]
    c = poly_add(poly_mul(neg_a, s), poly_add(m, e))
    return c, a


def decrypt(s, ct):
    # 相位 = c0 + c1·s（跟 4-1 的 c0 + <a, s> 一樣，只是換成多項式乘法），
    # 每個係數中心化，再 decode(中心化後的係數 list, DELTA)
    m = poly_centered(poly_add(ct[0], poly_mul(ct[1], s)))
    return decode(m, DELTA)


if __name__ == "__main__":
    s = keygen()
    z1 = [1.5, -2.0, 3.25, 0.5]
    z2 = [2.0, 0.5, -1.0, 4.0]

    enc1 = encrypt(s, z1)
    enc2 = encrypt(s, z2)

    # 1. 加解密
    print(f"z1={z1}, decrypt={decrypt(s, enc1)}")
    print(f"z2={z2}, decrypt={decrypt(s, enc2)}")

    # 2. 密文相加
    enc_add = tuple(poly_add(i, j) for i, j in zip(enc1, enc2))
    print(f"z1+z2={add(z1, z2)}, decrypt={decrypt(s, enc_add)}")

    # 3. 誤差比較：重複加解密 1000 次，記錄每個 slot 誤差的平均絕對值
    errors = []
    trails = 1000
    for _ in range(trails):
        z = [random.randint(-1000, 1000) / 100. for _ in range(N // 2)]
        errors += [abs(i - j) for i, j in zip(z, decrypt(s, encrypt(s, z)))]
    print(f"每個數字 平均誤差={np.mean(errors):.2e} 標準差={np.std(errors):.2e}")

    # 4.（選做）sample extraction
    z = [random.randint(-1000, 1000) / 100. for _ in range(N // 2)]
    c = encrypt(s, z)
    expected = sum(z)

    a_vec = [c[1][0]] + [-c[1][i] for i in range(N - 1, 0, -1)]
    m = (c[0][0] + sum(mul(a_vec, s)) + Q // 2) % Q - Q // 2
    result = m * N / 2 / DELTA
    print(f"Sample Extraction 預期 {expected:.2f} 結果 {result:.2f} 誤差 {abs(expected - result):.2e}")
