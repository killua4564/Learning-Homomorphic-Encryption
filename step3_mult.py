"""
步驟 3-2：密文相乘 —— 雜訊預算與乘法深度

目的：
  在 LWE 上做出「不需要私鑰」的密文乘法，觀察每次乘法讓雜訊長多少、
  q 的大小決定能連乘幾次（乘法深度），以及乘法讓密文長度爆炸的問題。

這裡換成「訊息放在低位」的編碼（BGV 風格），而且用對稱金鑰（只有私鑰，沒有公鑰），
讓乘法最容易看懂。這個實驗只看雜訊，不追求安全性，所以 n 取很小。

方案：
  參數：維度 n、明文模數 t=2、密文模數 q

  KeyGen：
    s  <- {-1, 0, 1}^n
    S  =  (1, s_1, ..., s_n)                 (解密用的「鑰匙向量」，長度 n+1)

  Enc(m)，m ∈ {0, 1}：
    a  <- Z_q^n 均勻隨機
    e  <- {-3, ..., 3}
    c  =  ( -<a, s> + m + t·e  mod q ,  a_1, ..., a_n )
    → 檢查：<c, S> = m + t·e  (mod q)，這個值叫「相位」(phase)

  Dec(c, S)：
    p  =  <c, S> mod q，再中心化到 (-q/2, q/2]
    m  =  p mod t
    → 只要 |m + t·e| < q/2，就能解對

  Add(c, c')：
    c + c'（逐項相加 mod q），鑰匙不變

  Mul(c, c')：outer product（張量積）
    <c, S> · <c', S'> = Σ_i Σ_j (c_i c'_j)(S_i S'_j) = <c ⊗ c', S ⊗ S'>
    → c'' = c ⊗ c' 只用到密文，做運算的一方不需要私鑰；
      S'' = S ⊗ S' 由持有私鑰的一方在解密時自己算
    → 代價：密文長度每乘一次就乘上 (n+1)，真實方案用 relinearization 壓回原長度

實驗：
  從 c = Enc(1) 開始，連續做 c = Mul(c, Enc(1))，直到解密出錯或達到 MAX_MULS 次；
  每個 q 跑 50 次，統計每次乘法增加的雜訊 bit 數，以及平均第幾次乘法出錯。

結果（n=2，一次執行的結果）：
  q = 2^8+1    平均第 4.3 次出錯
  q = 2^14+1   平均第 8.3 次出錯
  q = 2^20+1   平均第 10.6 次出錯（預測約 12.5 次，多數試驗跑到上限 11 次仍未出錯，故被低估）

觀察：
  - 每乘一個新密文，相位就乘上 |1 + 2e|，雜訊平均增加 E[log2|1 + 2e|] ≈ 1.52 bit
    （程式用 bit_length 取整數量測，數字會偏低）
  - 雜訊超過 q/2 就解錯 → log2(q/2) 就是「雜訊預算」，能撐的乘法次數約為 (log2 q − 1) / 1.52
  - q 必須與 t 互質：若 t | q（例如 q = 2^20），把密文 mod t 會讓雜訊 t·e 消失，
    LWE 退化成 F_2 上無雜訊的線性方程組，已知少量明文即可用高斯消去解出 s mod 2

執行：sage -python step3_mult.py
"""
import random

import numpy as np
from sage.all import Zmod, vector

n, t = 2, 2
MAX_MULS = 11  # 最多連乘幾次；密文長度為 (n+1)^(次數+1)，n=2 時 11 次約 53 萬


def keygen():
    s = F([random.randint(-1, 1) for _ in range(n)])
    S = vector([1] + s.list(), F.base_ring())
    return s, S


def encrypt(s, m):
    a = F.random_element()
    e = random.randint(-3, 3)
    return vector([-a * s + m + e * t] + a.list(), F.base_ring())


def decrypt(c, S):
    p = (int(c * S) + q // 2) % q - q // 2
    return p % t, abs(p).bit_length()


def add(c1, c2):
    return c1 + c2


def mul(c1, c2, S1, S2):
    c = vector(c1.outer_product(c2)).change_ring(F.base_ring())
    S = vector(S1.outer_product(S2)).change_ring(F.base_ring())
    return c, S


if __name__ == "__main__":
    for q in (2**8 + 1, 2**14 + 1, 2**20 + 1):  # q 必須與 t 互質；改成 2 的次方會出事（見上方說明）
        F = Zmod(q) ** n
        pk, sk = keygen()

        bits = []
        error_counts = []
        for _ in range(50):
            c, csk = encrypt(pk, 1), sk
            _, prev_bit_length = decrypt(c, csk)
            for count in range(1, MAX_MULS + 1):
                c, csk = mul(c, encrypt(pk, 1), csk, sk)
                bit, current_bit_length = decrypt(c, csk)
                bits.append(current_bit_length - prev_bit_length)
                prev_bit_length = current_bit_length
                if bit != 1:
                    error_counts.append(count)
                    break
            else:
                # 乘到上限都沒出錯：記成上限值（會讓平均略為低估）
                error_counts.append(MAX_MULS)

        print(f"{q=} 平均增加 {np.mean(bits):.2f} 個位元，平均第 {np.mean(error_counts):.2f} 次乘法出現錯誤")

