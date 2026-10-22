"""
步驟 2：LWE 加密（Regev 2005）與 LLL 攻擊

目的：
  (1) 實作最基本的 LWE 公鑰加密（一次加密 1 個位元），觀察解密雜訊與失敗率
  (2) 只用公鑰，以格基約化（LLL）回推私鑰，體會參數變大後攻擊為什麼失效

方案：
  參數：維度 n、模數 q、樣本數 m、雜訊範圍 B（e 的每一項從 {-B, ..., B} 均勻取）

  KeyGen：
    s  <- Z_q^n                       (私鑰)
    A  <- Z_q^{m x n} 均勻隨機
    e  <- {-B..B}^m
    b  =  A·s + e  (mod q)            (公鑰 = (A, b))

  Enc(bit μ)：
    x  <- {0,1}^m                     (隨機挑一個子集合)
    u  =  xᵀ·A        (mod q)         (長度 n 的向量)
    v  =  xᵀ·b + μ·⌊q/2⌋  (mod q)     (一個數)
    密文 = (u, v)

  Dec(u, v)：
    d  =  v − <u, s>  =  μ·⌊q/2⌋ + x·e  (mod q)
    d 比較靠近 0 就輸出 0，比較靠近 q/2 就輸出 1

  解密雜訊 x·e：最壞情況可達 m·B = 40 > q/4 ≈ 24，不保證永遠解對；
  但標準差只有 √(40/3) ≈ 3.65，q/4 約在 6.6 個標準差外，失敗機率約 10⁻¹⁰。

攻擊（Kannan embedding + LLL）：
  b − A·s ≡ e (mod q)，所以 (e, 1) 是下面這個格裡的短向量：
      [ q·I_m   0 ]
      [ Aᵀ      0 ]
      [ bᵀ      1 ]
  格的維度 d = m+1，行列式 det = q^(m−n)。
  LLL 找到 (e, ±1) 之後，解 A·s = b − e 即可得到 s。

  觀察（每組 5 次）：m=40 時 n ≤ 26 大多能解、n ≥ 28 失敗；m=80 可撐到 n ≈ 34。
  解釋：目標向量長度約 √(2m/3 + 1)，格的「正常最短向量」約為 Gaussian heuristic
        √(d/2πe)·det^(1/d)；兩者比值（間隙比）越大越容易抓到，而 LLL 能分辨的
        間隙約隨 δ^d 指數成長（δ ≈ 1.02，見 step2_lll_delta.py）。
        → n、B 變大會變難；q 變大會變簡單；m 存在對攻擊者最有利的值。

執行：sage -python step2_lwe.py
"""
import random

from sage.all import (
    GF,
    ZZ,
    Matrix,
    block_matrix,
    identity_matrix,
    random_matrix,
    random_vector,
    vector,
    zero_matrix,
)

n, q, m, B = 26, 97, 40, 1
F = GF(q)


def keygen():
    s = random_vector(F, n)
    A = random_matrix(F, m, n)
    e = vector([random.randint(-B, B) for _ in range(m)], F)
    b = A * s + e
    return (A, b), (s,)


def encrypt(pk, bit):
    (A, b) = pk
    x = vector(random.choices(range(2), k=m), F)
    u = x * A
    v = x * b + bit * q // 2
    return (u, v)


def decrypt(sk, ct):
    (s,), (u, v) = sk, ct
    d = v - u * s
    return round(int(d) / (q // 2)) & 1


if __name__ == "__main__":
    pk, sk = keygen()
    trials, fails = 10000, 0
    for _ in range(trials):
        bit = random.randrange(2)
        if decrypt(sk, encrypt(pk, bit)) != bit:
            fails += 1
    print(f"[+] {trials} 次加解密，失敗 {fails} 次")

    A, b = pk
    sol_matrix = block_matrix(ZZ, [
        [q * identity_matrix(m), zero_matrix(ZZ, m, 1)],
        [       Matrix(ZZ, A.T), zero_matrix(ZZ, n, 1)],
        [         Matrix(ZZ, b),     Matrix(ZZ, [[1]])],
    ])
    # print(sol_matrix)
    for sol in sol_matrix.LLL():
        if all(ei in (-1, 0, 1) for ei in sol[:-1]) and sol[-1] in (-1, 1):
            # print(sol)
            e = sol[-1] * vector(F, sol[:-1])
            s = A.solve_right(b - e)
            assert s == sk[0]
            print("[+] 私鑰回推成功")
            break
    else:
        print("[-] LLL 沒有找到 (e, ±1)，這次私鑰回推失敗（n 越接近 m 越容易失敗）")
