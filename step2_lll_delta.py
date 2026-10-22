"""
步驟 2 補充：自己量出 LLL 的「1.0219」

做法跟 Gama–Nguyen (EUROCRYPT 2008) 一樣：
  1. 產生隨機格（Goldstein–Mayer 形式：第一列 (p, 0, ..., 0)，其餘列 (x_i, 單位向量)）
  2. 跑 LLL，拿第一個（最短的）向量 b1
  3. 算 root Hermite factor：δ = (‖b1‖ / det^(1/d))^(1/d)
     也就是假設 ‖b1‖ ≈ δ^d · det^(1/d)，反解出 δ

執行：sage -python step2_lll_delta.py   （d=160 要跑大約 40 秒）
"""
from sage.all import ZZ, RR, Matrix, block_matrix, exp, identity_matrix, log, random_prime, randrange

for d in (40, 80, 120, 160):
    # det = p，取 p 大約 10·d bits，讓格「夠隨機」
    p = random_prime(2 ** (10 * d), lbound=2 ** (10 * d - 1))
    col = Matrix(ZZ, d - 1, 1, [randrange(p) for _ in range(d - 1)])
    M = block_matrix(ZZ, [
        [Matrix(ZZ, [[p]]), Matrix(ZZ, 1, d - 1)],
        [col, identity_matrix(ZZ, d - 1)],
    ])
    b1 = M.LLL()[0]
    # 用 log 計算，避免 det^(1/d) 這種超大數字溢位
    delta = exp((log(RR(b1.norm())) - log(RR(p)) / d) / d)
    print(f"d={d:4d}  δ ≈ {delta:.5f}")
