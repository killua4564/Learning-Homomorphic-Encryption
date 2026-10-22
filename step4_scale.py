"""
步驟 4-1：把實數放進 LWE —— scale Δ 與精度

目的：
  LWE 的相位只能是整數。CKKS 的作法是把實數 x 乘上 scale Δ、四捨五入成整數再加密，
  雜訊直接加在這個整數的最低位，解密後除以 Δ，得到「近似」的 x。
  這支程式量測 Δ 與精度的關係，以及 Δ 的上限。

方案（步驟 3-2 的對稱金鑰版本，但沒有明文模數 t）：
  KeyGen：s ← {−1, 0, 1}^n
  Enc(x)：a ← Z_q^n 均勻隨機，e ← {−3, …, 3}
          c = ( −<a, s> + round(Δ·x) + e  mod q ,  a_1, …, a_n )
          → 相位 = c_0 + <a, s> = round(Δ·x) + e   (mod q)
  Dec(c)：p = 相位中心化到 (−q/2, q/2]；x' = p / Δ

實驗：n = 2、q ≈ 2^40、x = π，每個 Δ 加解密 1000 次，記錄 |x' − x| 的平均與最大值。

結果（一次執行）：
  Δ        2^4      2^8      2^16     2^24     2^36      2^37      2^38   2^40
  平均誤差 1.1e-1   6.8e-3   2.6e-5   1.1e-7   2.6e-11   1.3e-11   4.00   3.00

觀察：
  - 誤差 ≈ 雜訊 / Δ：Δ 每大 2^8 倍，誤差就小 2^8 倍。雜訊大小由安全性決定，
    使用者只能透過 Δ（TenSEAL 的 global_scale）選擇精度。
  - 必須 Δ·|x| < q/2，即 log Δ < log q − 1 − log|x|（這裡約 2^37.35）。
    超過後相位會繞回：解出 x − k·(q/Δ)，k 是讓相位落回 (−q/2, q/2] 所需扣掉的 q 的次數。
      Δ = 2^38：q/Δ ≈ 4，k = 1 → x' ≈ −0.86，誤差 4.00
      Δ = 2^40：q/Δ ≈ 1，k = 3 → x' ≈ 0.14， 誤差 3.00
  - 因此 CKKS 對雜訊是「精度逐漸變差」，對溢位則是「懸崖式出錯」。
  - 對應到真實參數：coeff_mod_bit_sizes 的第一個質數要裝得下 Δ（精度）與 |x|（整數部分）。

執行：sage -python step4_scale.py
"""
import random

import numpy as np
from sage.all import (
    GF,
    VectorSpace,
    vector,
)

n = 2
q = 2**40 + 15          # 約 2^40 的奇數（這一步不需要質數）
x = 3.14159265358979    # 要加密的實數

V = VectorSpace(GF(q), n)


def keygen():
    # 回傳 s（長度 n 的 list，每項取自 {−1, 0, 1}）
    return vector([random.randint(-1, 1) for _ in range(n)], V.base_ring())


def encrypt(s, x, delta):
    # 回傳 c = [c_0, a_1, ..., a_n]
    a = V.random_element()
    e = random.randint(-3, 3)
    return vector([-a * s + round(delta * x) + e] + a.list(), V.base_ring())


def decrypt(s, c, delta):
    # 回傳 x'（浮點數）
    S = vector([1] + s.list(), V.base_ring())
    p = (int(c * S) + q // 2) % q - q // 2
    return p / delta


if __name__ == "__main__":
    s = keygen()
    trails = 1000

    # 1. 先跑一次看看
    delta = 2**16
    x_dec = decrypt(s, encrypt(s, x, delta), delta)
    print(f"Δ=2^16：x' = {x_dec}，誤差 = {abs(x_dec - x):.2e}")

    # 2. 改變 Δ，每組跑 1000 次
    for k in (4, 8, 16, 24):
        delta = 2**k
        errors = [
            abs(decrypt(s, encrypt(s, x, delta), delta) - x)
            for _ in range(trails)
        ]
        print(f"Δ=2^{k:<2}  平均誤差 {np.mean(errors):.2e}  最大誤差 {max(errors):.2e}")

    # 3. 試 Δ = 2^36、2^37、2^38、2^40，觀察結果
    for k in (36, 37, 38, 40):
        delta = 2**k
        errors = [
            abs(decrypt(s, encrypt(s, x, delta), delta) - x)
            for _ in range(trails)
        ]
        print(f"Δ=2^{k:<2}  平均誤差 {np.mean(errors):.2e}  最大誤差 {max(errors):.2e}")
