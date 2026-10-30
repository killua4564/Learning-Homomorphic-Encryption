"""
步驟 5-0a：CKKS 編碼 —— 把一串數字變成一個整數係數多項式（還不加密）

目的：
  TenSEAL 的 ckks_vector(ctx, [...]) 會把一整串數字裝進「一個」多項式再加密。
  這一步只做「裝進多項式」與「拿出來」，看懂為什麼：
      多項式相加 → 每個 slot 各自相加
      多項式相乘（mod X^N + 1）→ 每個 slot 各自相乘

方法：
  環 Z[X]/(X^N + 1)，N = 8。X^N + 1 的 N 個複數根是 ζ^1, ζ^3, ..., ζ^15，ζ = e^(iπ/N)。
  多項式 m(X) 有 N 個實數係數；它在 N 個根上的值兩兩共軛，所以只能自由指定 N/2 = 4 個值。
  這 4 個值就是 4 個 slot。這裡選的 4 個點是 ζ^(5^j mod 2N)（CKKS 原始論文的選法；SEAL 用 3 的次方），j = 0..3：
      ζ^1, ζ^5, ζ^9, ζ^13       （另外 4 個根 ζ^15, ζ^11, ζ^7, ζ^3 是它們的共軛）
  （為什麼是 5 的次方，5-0e 做旋轉時會用到。）

  encode(z, Δ)：
    1. 把 4 個 slot 值 z 和它們的共軛排成 N 個「目標值」，對應到 N 個根
    2. 解線性方程組 V · m = 目標值（V[i][k] = 根_i 的 k 次方），得到 N 個係數
    3. 係數乘上 Δ、四捨五入 → 整數係數多項式（這才是之後要加密的東西）
  decode(m, Δ)：把多項式代入 4 個 slot 點，除以 Δ

結果：
  z1 = [1.5, −2, 3.25, 0.5] 編碼成 8 個整數係數 [851968, −173036, 1158524, −781030, 0, 781030, −1158524, 173036]
  - decode(encode(z)) ≈ z，誤差約 1e-6（四捨五入到 Δ = 2^20 造成）
  - 多項式相加 → slot 各自相加；多項式相乘（mod X^N + 1）→ slot 各自相乘（decode 時除以 Δ²）
  - 係數後半段是前半段的「負號 + 倒序」，虛部全為 0：共軛對稱讓 N 個實數係數剛好對應 N/2 個 slot

觀察：
  - slot 點都是 X^N + 1 的根（ζ^N = −1），在這些點上把 X^N 換成 −1 不改變值，
    所以只有乘法需要 mod X^N + 1，加法和求值都不受影響
  - sum(z1) = 3.25 = 4 × 0.8125 = (N/2)·m[0]/Δ：所有根的 k 次方加總在 k ≠ 0 時為 0，
    只剩常數項。推廣成 m_k = (1/N)·Σ_j target_j·conj(root_j)^k，就是複數版的 intt
  - 每個 slot 其實是複數；把兩個實數塞進同一個 slot 的實部和虛部，加法沒問題，
    但乘法會變成複數乘法、實部混進虛部的乘積，所以 TenSEAL 只用實部

執行：python3 step5_0a_encode.py
"""
import numpy as np

N = 8
DELTA = 2**20
ZETA = np.exp(1j * np.pi / N)  # 2N 次本原單位根

# 4 個 slot 點：[ZETA**1, ZETA**5, ZETA**9, ZETA**13]
SLOT_ROOTS = np.array([ZETA ** pow(5, j, 2 * N) for j in range(N // 2)])
SLOT_ROOTS_MATRIX = np.array([SLOT_ROOTS ** i for i in range(N)]).transpose()

# 全部 N 個根：先放 4 個 slot 點，再放它們的共軛（順序對應好，encode 時比較方便）
ALL_ROOTS = np.concatenate([SLOT_ROOTS, np.conj(SLOT_ROOTS)])
ALL_ROOTS_MATRIX = np.array([ALL_ROOTS ** i for i in range(N)]).transpose()


def add(a, b):
    return [i + j for i, j in zip(a, b)]


def mul(a, b):
    return [i * j for i, j in zip(a, b)]


def negacyclic_mul(a, b):
    # 整數係數多項式相乘 mod (X^N + 1)，直接展開（步驟 3-3 的 schoolbook 版本，不 mod p）
    c = [0] * N
    for i in range(N):
        for j in range(N):
            if i + j < N:
                c[i + j] += a[i] * b[j]
            else:
                c[i + j - N] -= a[i] * b[j]
    return c


def encode(z, delta):
    # target = z 接上 z 的共軛（長度 N），對應 ALL_ROOTS
    # V[i][k] = ALL_ROOTS[i] ** k
    # m = np.linalg.solve(V, target)  → 虛部應該 ≈ 0
    # 回傳 [round(delta * 實部) for ...]，用 Python int 的 list
    target = np.concatenate([z, np.conj(z)])
    m = np.linalg.solve(ALL_ROOTS_MATRIX, target)
    return [int(v) for v in (m.real * delta).round()]


def decode(m, delta):
    # 把多項式 m 代入 SLOT_ROOTS 的每一個點，除以 delta，回傳 4 個複數
    m = np.array(m, dtype=np.float64).reshape(N, 1)
    z = np.matmul(SLOT_ROOTS_MATRIX, m)
    return [float(v) for v in (z.real / delta).reshape(N // 2)]


if __name__ == "__main__":
    z1 = [1.5, -2.0, 3.25, 0.5]
    z2 = [2.0, 0.5, -1.0, 4.0]

    m1, m2 = encode(z1, DELTA), encode(z2, DELTA)
    print("m1 的 8 個整數係數：", m1)

    # 1. 來回一次
    # 印出 decode(m1) 的實部，和 z1 比較
    print(f"z1={z1}, decode={decode(m1, DELTA)}")

    # 2. 相加
    # 逐項相加 m1 + m2，decode，和 z1 + z2 比較
    print(f"z1+z2={add(z1, z2)}, decode={decode(add(m1, m2), DELTA)}")

    # 3. 相乘
    # negacyclic_mul(m1, m2)，decode 時除以 DELTA**2，和 z1 * z2 比較
    print(f"z1*z2={mul(z1, z2)}, decode={decode(negacyclic_mul(m1, m2), DELTA ** 2)}")

    # 4. slot 總和 vs 常數項
    # 印出 sum(z1) 和 m1[0] / DELTA，找出兩者的關係
    print(f"sum(z1)={sum(z1)}, m1[0]/DELTA={m1[0] / DELTA}")
