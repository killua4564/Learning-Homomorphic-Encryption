"""
步驟 3-1：LWE 密文連加的雜訊成長

目的：
  觀察同態「加法」會讓雜訊怎麼長大、加到第幾個密文時解密開始出錯，
  並用數學推導預測雜訊的標準差，再用實驗驗證。

方案：
  沿用步驟 2 的 Regev LWE 位元加密（規格見 step2_lwe.py），參數 n=26, q=97, m=40, B=1。
    密文 (u, v)，解密時 d = v − <u, s> = μ·⌊q/2⌋ + x·e  (mod q)
  同態加法：(u1 + u2, v1 + v2)，解密結果是兩個位元的 XOR（mod 2 的加法）。

實驗：
  對 k = 1, 2, 5, 10, 20, 40，各做 10000 次試驗：
    1. 產生一把新金鑰
    2. 隨機產生 k 個位元，加密後把 k 個密文全部加起來
    3. 解密，與 k 個位元的 XOR 比對，記錄失敗率
    4. 記錄雜訊 = d − (XOR 結果)·⌊q/2⌋，中心化到 [−48, 48]，統計平均值與標準差

  開關 ZERO_COV 決定加密時隨機向量 x 的取法：
    False：x ← {0, 1}^m        （Regev 原始版本，E[x_i] = 1/2）
    True ：x ← {−1, 0, 1}^m    （以 0 為中心，E[x_i] = 0）

預測公式：
  標準差
    ZERO_COV=False：Var = k·(40/3) + 2·C(k,2)·(40/6) = (20/3)·k·(k+1)
        → 同一把金鑰的 k 個密文共用同一組 e，彼此有共變異數，雜訊約跟 k 成長
    ZERO_COV=True ：Var = k·m·E[x²]·E[e²] = k·40·(2/3)·(2/3) = k·160/9
        → E[x_i] = 0 讓共變異數消失，雜訊回到跟 √k 成長
  平均值
    q=97 是奇數，2·⌊q/2⌋ = 96 ≡ −1 (mod q)，每兩個 1 相加會偏移 −1，
    k 個隨機位元的平均偏移為 −(k − 1)/4

結果摘要（完整數據與分析見 README.md「加法雜訊」一節）：
  - 兩種 x 的實測標準差都與預測吻合，直到雜訊接近飽和（約 28，即 [−48, 48] 均勻分布的標準差）
  - k=10 時失敗率：ZERO_COV=False 約 36%，ZERO_COV=True 約 7.5%
  - 雜訊飽和後失敗率約 50%：密文已不含任何明文資訊

執行：sage -python step3_noise.py
"""
import functools
import math
import operator
import random

import numpy as np
from sage.all import GF, random_matrix, random_vector, vector

n, q, m, B = 26, 97, 40, 1
F = GF(q)

ZERO_COV = True


def keygen():
    s = random_vector(F, n)
    A = random_matrix(F, m, n)
    e = vector([random.randint(-B, B) for _ in range(m)], F)
    b = A * s + e
    return (A, b), (s,)


def encrypt(pk, bit):
    (A, b) = pk
    if ZERO_COV:
        x = vector([random.randint(-1, 1) for _ in range(m)], F)
    else:
        x = vector(random.choices(range(2), k=m), F)
    u = x * A
    v = x * b + bit * q // 2
    return (u, v)


def decrypt(sk, ct):
    (s,), (u, v) = sk, ct
    d = v - u * s
    return round(int(d) / (q // 2)) & 1, int(d)


def add_cipher(c1, c2):
    (u1, v1), (u2, v2) = c1, c2
    return (u1 + u2), (v1 + v2)


if __name__ == "__main__":
    for k in (1, 2, 5, 10, 20, 40):
        trials, fails, errors = 10000, 0, []
        for _ in range(trials):
            pk, sk = keygen()
            bits = [random.randrange(2) for _ in range(k)]
            c = functools.reduce(add_cipher, [encrypt(pk, bit) for bit in bits])

            mbit, d = decrypt(sk, c)
            bit = functools.reduce(operator.xor, bits)
            fails += mbit != bit

            error = d - bit * q // 2
            errors.append((error + q // 2) % q - q // 2)

        print(f"[+] {k=} {ZERO_COV=}")
        print(f"{trials} 次加解密，失敗率 {100. * fails / trials:.2f}%")
        if ZERO_COV:
            print(f"預測：平均值 {-(k - 1) / 4:.2f}、標準差 {math.sqrt(k * 160 / 9):.2f}")
        else:
            print(f"預測：平均值 {-(k - 1) / 4:.2f}、標準差 {math.sqrt(20 * k * (k + 1) / 3):.2f}")
        print(f"實際：平均值 {np.mean(errors):.2f}、標準差 {np.std(errors):.2f}")
