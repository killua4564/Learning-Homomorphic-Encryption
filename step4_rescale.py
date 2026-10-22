"""
步驟 4-3：rescale 與模數鏈 —— 把倍數拉回 Δ

目的：
  實作 CKKS 的 rescale，觀察它對乘法次數與「不同深度的密文能否相加」的影響，
  最後示範玩具參數下只用一個密文就能被 LLL 找回私鑰。

方案：
  rescale = 把整個密文除以一個質數 p 並四捨五入：
      相位 <c, S> = P + k·q
      → <round(c / p), S> ≈ P/p + k·(q/p)      （外加一點四捨五入誤差）
  倍數 Δ² 變成 Δ²/p ≈ Δ，模數從 q 變成 q/p。為了讓 q/p 是整數，q 取成質數的乘積（模數鏈）：

      PRIMES = [ base,     p1,      p2,      p3    ]
                 約 2^25   約 2^20  約 2^20  約 2^20
      密文一開始的模數 q = base·p1·p2·p3；每 rescale 一次就從最後面拿掉一個質數。
      base 留到最後，用來裝 Δ·x（Δ 的 20 bit + x 的整數部分 + 餘裕）。

  密文表示為 (c, scale, level)，level = 目前還剩幾個質數，鑰匙 S 另外傳遞。
  - add / mul：兩邊 level 不同時，取較低的 level，直接 mod 較小的模數對齊
    （較小的模數整除較大的模數，所以相位不變）
  - lift：把鑰匙是 S 的密文補 0，改寫成鑰匙是 S⊗S 的密文，才能與乘過的密文相加
    （S⊗S 的前 len(S) 項正好是 S；真實方案改用 relinearization 把 S⊗S 壓回 S）

  ※ 真實方案的質數還要 ≡ 1 (mod 2N)（NTT 需要，見步驟 3-3），並多一個只在
    relinearization 使用的特殊質數；這個玩具都省略了。

實驗與結果：
  A. x = 1.5 連續「乘上新加密的 1.5 再 rescale」
     → 前 3 次相對誤差約 1e-6，scale 一直維持約 2^20；第 4 次 base 也被用掉（level 0），解出 0。
     → 能乘的次數 = 中間質數的個數（3）。與總 bit 數相同（85 bit）、不 rescale 的版本相比，
       次數一樣：rescale 不會多出乘法空間。
  B. 計算 x² + x（正確值 3.75）
     (1) x² 不 rescale 直接加上 x：誤差 1.50 —— x 的倍數比 x² 小一個 Δ，解密除以 Δ² 後 x 消失
     (2) x² 先 rescale 再加上 x：
         質數離 2^20 約 0.2%（註解掉的那組 PRIMES）→ 誤差約 2.9e-3
         質數離 2^20 只差 3～17（目前這組）       → 誤差約 2e-5
     → rescale 真正的用途是「對齊倍數」；剩下的誤差來自 p ≠ Δ 造成的 scale 不一致，
       質數越接近 Δ，誤差越小，最後被加密雜訊蓋過。SEAL 也因此挑非常接近 2^bit 數的質數，
       並在 scale 不一致時拒絕相加（scale mismatch）。
  C. 只用一個密文 (c0, c1, c2) 做 LLL 攻擊（Kannan embedding）
         [ c0  W  0  0 ]
         [ c1  0  W  0 ]
         [ c2  0  0  W ]
         [ q   0  0  0 ]
     目標向量 ±(Δx + e, W, s1·W, s2·W)，取 W = 2^21 讓各座標大小相近。
     格的 det = q·W^3 ≈ 2^148，Gaussian heuristic 約 2^36，目標約 2^22，間隙極大，
     LLL 每次都能找回 s。真實參數中一個 RLWE 密文相當於 N 組 LWE 樣本，維度數千，無法這樣攻擊。

執行：sage -python step4_rescale.py
"""
import math
import random

from sage.all import (
    block_matrix,
    identity_matrix,
    vector,
    zero_matrix,
)

n = 2
x = 1.5
DELTA = 2**20
# PRIMES = [33554393, 1048573, 1047559, 1046557]  # 對照組：離 2^20 約差 0.2%，B(2) 誤差約 2.9e-3
PRIMES = [33554393, 1048573, 1048571, 1048559]   # 離 2^20 分別只差 3、5、17

def modulus(level):
    # level = 剩下幾個質數；模數 = 前 level 個質數的乘積
    q = 1
    for p in PRIMES[:level]:
        q *= p
    return q


def keygen():
    return vector([random.randint(-1, 1) for _ in range(n)])


def encrypt(s, x):
    # 跟 4-2 一樣，模數用完整的模數鏈；回傳 (c, scale, level)
    q = modulus(len(PRIMES))
    a = vector([random.randint(0, q-1) for _ in range(n)])
    e = random.randint(-3, 3)
    c = vector([(-a * s + round(DELTA * x) + e) % q] + a.list())
    return c, DELTA, len(PRIMES)


def decrypt(ct, S):
    # 用這個密文自己 level 的模數中心化相位，再除以 scale
    c, scale, level = ct
    q = modulus(level)
    p = (int(c * S) + q // 2) % q - q // 2
    return p / scale


def lift(ct, S):
    # 把鑰匙是 S 的密文，改寫成鑰匙是 S⊗S 的密文（相加時兩邊鑰匙要一樣）
    # 原理：S⊗S 的前 len(S) 項是 1·S_j = S_j，所以 c 後面補 0 即可，相位不變
    # （真實方案用 relinearization 反過來把 S⊗S 壓回 S）
    c, scale, level = ct
    c = vector(c.list() + [0] * (len(S) * len(S) - len(S)))
    S = vector(S.outer_product(S))
    return (c, scale, level), S


def add(ct1, S1, ct2, S2):
    (c1, scale1, level1), (c2, scale2, level2) = ct1, ct2
    assert S1 == S2, "相加的兩個密文鑰匙要一樣"
    level = min(level1, level2)
    q = modulus(level)
    return ((c1 + c2) % q, scale1, level), S1  # 注意：scale 直接沿用 ct1 的


def mul(ct1, S1, ct2, S2):
    # 兩邊取較低的 level（mod 較小的模數），再做 outer product；scale 相乘
    (c1, scale1, level1), (c2, scale2, level2) = ct1, ct2
    level = min(level1, level2)
    q = modulus(level)
    c = vector((c1 % q).outer_product(c2 % q))
    S = vector(S1.outer_product(S2))
    return (c, scale1 * scale2, level), S


def rescale(ct, S):
    # p = 目前 level 的最後一個質數
    # c 的每一項：先中心化，再除以 p 四捨五入，再 mod 新的模數
    # scale 除以 p，level 減 1
    c, scale, level = ct

    p = PRIMES[level - 1]
    q = modulus(level)
    c = vector([(((i + q // 2) % q - q // 2) / p).round() for i in c])
    c = c % (q // p)

    return (c, scale / p, level - 1), S


if __name__ == "__main__":
    s = keygen()

    print("A. 連續相乘 + rescale")
    ct, expected = encrypt(s, x), x
    skey = vector([1] + s.list())
    result = decrypt(ct, skey)
    print(f"第 0 次乘法，誤差 {abs(result - expected):.2e}，log2(scale)={math.log2(ct[1]):.2f}，log2(expected)={math.log2(expected):.2f}")
    for i in range(1, len(PRIMES) + 1):
        # 乘上新加密的 x 再 rescale；最後一次會把 base 也用掉（level 0），解出 0
        ct, skey = rescale(*mul(ct, skey, encrypt(s, x), vector([1] + s.list())))
        expected = expected * x
        result = decrypt(ct, skey)
        print(f"第 {i} 次乘法，誤差 {abs(result - expected):.2e}，log2(scale)={math.log2(ct[1]):.2f}，log2(expected)={math.log2(expected):.2f}")

    print("\nB. x² + x")
    cx = encrypt(s, x)
    skey = vector([1] + s.list())
    cx2, skey2 = mul(cx, skey, encrypt(s, x), skey)
    # (1) x² 不 rescale，直接加上 x
    r1 = decrypt(*add(cx2, skey2, *lift(cx, skey)))
    print(f"1. 誤差 {abs(r1 - x ** 2 - x):.2e}")
    # (2) x² 先 rescale，再加上 x（add 會自動把 x 降到同一個 level）
    r2 = decrypt(*add(*rescale(cx2, skey2), *lift(cx, skey)))
    print(f"2. 誤差 {abs(r2 - x ** 2 - x):.2e}")

    print("\nC. LLL 攻擊：只用一個密文 cx 和 q 找回私鑰")
    W = 2**21
    q = modulus(cx[2])
    B = block_matrix([
        [                 cx[0], W * identity_matrix(3)],
        [q * identity_matrix(1),      zero_matrix(1, 3)],
    ])
    for sol in B.LLL():
        if abs(sol[1]) == W and all(abs(i) <= W for i in sol[2:]):
            assert sol[1] * vector(sol[2:]) / W**2 == s
            print("[+] 成功找到金鑰")
