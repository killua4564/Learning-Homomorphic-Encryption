"""
步驟 1：部分同態加密 —— Paillier 的加法同態

目的：
  從 RSA 的乘法同態（Enc(a)·Enc(b) = Enc(a·b)，RSA LSB oracle 攻擊利用的性質）出發，
  實作 Paillier，驗證它的加法同態，並理解「部分同態」為什麼只能做一種運算。

方案（g = n + 1 的簡化版）：
  KeyGen：p, q 為 1024-bit 質數，n = p·q，λ = lcm(p−1, q−1)，μ = L(g^λ mod n²)⁻¹ mod n
          其中 L(x) = (x − 1) / n；公鑰 (n, g)，私鑰 (λ, μ)
  Enc(m)：c = g^m · r^n mod n²，r 為與 n 互質的隨機數
  Dec(c)：m = L(c^λ mod n²) · μ mod n
  關鍵：(1 + n)^m ≡ 1 + m·n (mod n²)，明文被放在「指數」上，又能被線性地讀出來

同態性質：
  加法：     Enc(m1) · Enc(m2) mod n² = Enc(m1 + m2)
  純量乘法： Enc(m)^k mod n²       = Enc(k·m)    （k 為公開整數）

觀察：
  - 明文放在指數上，所以「密文乘 ↔ 明文加」；要做到明文相乘，等於要從 g^a、g^b 算出 g^(ab)，
    這是 Diffie–Hellman（CDH）形狀的難題 → 只能做一種運算
  - RSA 和 Paillier 的安全性都建立在因數分解上，Shor 演算法（量子求週期）可以同時破解兩者

執行：python3 step1_paillier.py   （需要 pycryptodome）
"""
import math

from Crypto.Random.random import getrandbits
from Crypto.Util.number import getPrime


def keygen():
    p, q = getPrime(1024), getPrime(1024)
    n = p * q
    lamb = math.lcm(p - 1, q - 1)

    g = n + 1
    mu = pow((pow(g, lamb, n ** 2) - 1) // n, -1, n)

    return (n, g), (lamb, mu)


def encrypt(pk, m):
    n, g = pk
    assert 0 <= m < n
    r = getrandbits(1000)
    while not (0 < r < n and math.gcd(r, n) == 1):
        r = getrandbits(1000)
    return pow(g, m, n ** 2) * pow(r, n, n ** 2) % n ** 2


def decrypt(pk, sk, c):
    (n, _), (lamb, mu) = pk, sk
    return (mu * (pow(c, lamb, n ** 2) - 1) // n) % n


def add_cipher(pk, c1, c2):
    (n, _) = pk
    return c1 * c2 % n ** 2


if __name__ == "__main__":
    pk, sk = keygen()
    for _ in range(10):
        m1, m2 = getrandbits(1000), getrandbits(1000)
        c1, c2 = encrypt(pk, m1), encrypt(pk, m2)
        assert decrypt(pk, sk, add_cipher(pk, c1, c2)) == m1 + m2
    print("[+] 加法同態驗證通過")

    # 純量乘法：不解密，只用公開的整數 k，做出 Enc(k·m) = Enc(m)^k
    m, k = getrandbits(1000), getrandbits(10)
    assert 0 <= m * k < pk[0]
    c = encrypt(pk, m)
    assert decrypt(pk, sk, pow(c, k, pk[0] ** 2)) == k * m
    print("[+] 乘法同態驗證通過")
