"""
步驟 4-2：實數密文相乘（不 rescale）—— 倍數 Δ 會怎麼長大

目的：
  觀察兩個 CKKS 式密文相乘後，倍數從 Δ 變成 Δ²、Δ³…，以及在同一個 q 底下
  Δ 的大小如何同時影響精度與能乘的次數。

方案：
  沿用 4-1 的加密（相位 = Δ·x + e），乘法用步驟 3-2 的 outer product：
      (Δ·x1 + e1)(Δ·x2 + e2) = Δ²·x1·x2 + Δ(x1·e2 + x2·e1) + e1·e2
  所以解密時要除以 Δ²。每個密文因此要記住自己「目前的倍數」scale：
      剛加密 scale = Δ，乘 k 次後 scale = Δ^(k+1)
  （TenSEAL 的密文物件裡也存著 scale。）

  密文表示為 (c, scale)，鑰匙 S 另外傳遞：剛加密時 S = (1, s)，相乘後 S = S1 ⊗ S2。

實驗：
  q ≈ 2^60，從 x = 1.5 開始連續乘上新加密的 1.5；
  Δ = 2^5、2^10、2^20，相對誤差超過 50% 就視為失敗。

結果：
  Δ       每次乘法吃掉約   預測可乘次數   實測第幾次壞掉
  2^5     5.6 bit           9              10
  2^10    10.6 bit          4              5
  2^20    20.6 bit          1              2

觀察：
  - 條件是 log2(scale) + log2|值| < log2(q/2)：倍數和數值本身都會占用 q 的 bit。
    （把 x 改成 1.0 時數值不占 bit，Δ = 2^5 就能撐到 scale = 2^55。）
  - 同一個 q 底下，Δ 小 → 能乘很多次但很不準（2^5 時一開始就有數 % 誤差）；
    Δ 大 → 很準但只能乘幾次。
  - 溢位後解出的值 ≈ 0（相位最多 ±q/2，除以超過 q 的 scale 後趨近 0），
    所以誤差會剛好等於正確值本身；密文自己不知道已經壞了，TenSEAL 則會主動檢查並
    丟出 scale out of bounds。
  - 停止條件若用「相對誤差 > 100%」，會因為溢位後誤差黏在 100% 附近而時停時不停，
    所以改用 50%。

執行：sage -python step4_mul.py
"""
import math
import random

from sage.all import (
    GF,
    VectorSpace,
    vector,
)

n = 2
q = 2**60 + 33          # 約 2^60 的奇數
x = 1.5

V = VectorSpace(GF(q), n)
MAX_MULS = 11           # 密文長度 = 3^(次數+1)，11 次約 53 萬，還跑得動


def keygen():
    return vector([random.randint(-1, 1) for _ in range(n)], V.base_ring())


def encrypt(s, x, delta):
    # 跟 4-1 一樣，額外回傳目前的倍數 scale = delta
    a = V.random_element()
    e = random.randint(-3, 3)
    c = vector([-a * s + round(delta * x) + e] + a.list(), V.base_ring())
    return c, delta


def decrypt(ct, S):
    # 相位 = <c, S> 中心化後，除以這個密文自己的 scale
    c, scale = ct
    p = (int(c * S) + q // 2) % q - q // 2
    return p / scale


def mul(ct1, ct2, S1, S2):
    # c 和 S 各自做 outer product，倍數相乘
    (c1, scale1), (c2, scale2) = ct1, ct2
    c = vector(c1.outer_product(c2))
    S = vector(S1.outer_product(S2))
    return (c, scale1 * scale2), S


if __name__ == "__main__":
    s = keygen()
    for k in (5, 10, 20):
        delta = 2**k
        skey = vector([1] + s.list(), V.base_ring())
        ct, expected = encrypt(s, x, delta), x
        print(f"Δ=2^{k}")
        result = decrypt(ct, skey)
        print(f"第 0 次乘法，誤差 {abs(result - expected):.2e}，log2(scale)={math.log2(ct[1]):.2f}，log2(expected)={math.log2(expected):.2f}")
        for i in range(1, MAX_MULS + 1):
            # 乘上新加密的 x，印出誤差、log2(scale) 與 log2(正確值)；相對誤差 > 50% 視為失敗
            ct, skey = mul(ct, encrypt(s, x, delta), skey, vector([1] + s.list(), V.base_ring()))
            expected = expected * x
            result = decrypt(ct, skey)
            print(f"第 {i} 次乘法，誤差 {abs(result - expected):.2e}，log2(scale)={math.log2(ct[1]):.2f}，log2(expected)={math.log2(expected):.2f}")
            if 2 * abs(result - expected) > expected:
                break

