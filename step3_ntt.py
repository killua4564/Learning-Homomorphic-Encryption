"""
步驟 3-3：負循環 NTT —— 在 Z_p[X]/(X^N + 1) 上做 O(N log N) 的多項式乘法

目的：
  RLWE 把 LWE 的「數字」換成多項式環 Z_q[X]/(X^N + 1) 的元素；
  這個環上的乘法等於「負循環矩陣 × 向量」，直接算是 O(N²)，用 NTT 可降到 O(N log N)。

方法：
  - ntt(a, w)：遞迴 Cooley–Tukey，回傳 [A(w^0), A(w^1), ..., A(w^(n−1))]
      A(x) = E(x²) + x·O(x²)，再用 w^(n/2) = −1 做蝴蝶運算：
      A(w^i) = E_i + w^i·O_i，A(w^(i+n/2)) = E_i − w^i·O_i
  - X^N + 1 的根是 ψ 的奇數次方（ψ 為 2N 次本原單位根），所以需要 2N | p − 1
      作法：先把係數 a_i 乘上 ψ^i（twist），再用 w = ψ² 做一般 NTT，
            即等於在 ψ, ψ³, ..., ψ^(2N−1) 上求值
  - intt：用 ψ⁻² 做一般 NTT，第 i 項乘上 ψ^(−i)，最後乘 N⁻¹
  - negacyclic_mul：twist → NTT → 逐點相乘 → INTT，全程長度 N，不需補零到 2N

參數：p = 12289（NewHope 用過的質數，p − 1 = 2^12 · 3），11 是 F_p* 的原根，N = 256

結果（一次執行）：N = 256 時，NTT 約 0.55 毫秒，直接展開約 3.2 毫秒（約快 6 倍），
                  N 越大差距越大（O(N log N) vs O(N²)）

與同態加密的關係：這就是 SEAL／TenSEAL 要求 coeff_mod_bit_sizes 裡每個質數 q
                  都滿足 q ≡ 1 (mod 2N) 的原因

執行：python3 step3_ntt.py
"""
import random
import time

p = 12289
g = 11
N = 256

INV_N = pow(N, -1, p)


def schoolbook_negacyclic(a, b):
    # 直接展開 O(N^2)，遇到 X^N 就換成 -1：用來驗證 NTT 結果、並比較速度
    c = [0] * N
    for i in range(N):
        for j in range(N):
            k = i + j
            if k < N:
                c[k] = (c[k] + a[i] * b[j]) % p
            else:
                c[k - N] = (c[k - N] - a[i] * b[j]) % p
    return c


def ntt(a, w):
    n = len(a)
    if n == 1:
        return a

    w2 = pow(w, 2, p)
    evens, odds = ntt(a[::2], w2), ntt(a[1::2], w2)

    x = 1
    result = [0] * n
    for i in range(n // 2):
        result[i] = (evens[i] + x * odds[i]) % p
        result[i + n // 2] = (evens[i] - x * odds[i]) % p
        x = x * w % p

    return result


def intt(A, psi):
    ipsi = pow(psi, -1, p)
    ipsi2 = pow(ipsi, 2, p)

    x = 1
    result = []
    for i, a in enumerate(ntt(A, ipsi2)):
        result.append(INV_N * x * a % p)
        x = x * ipsi % p

    return result


def negacyclic_mul(a, b):
    psi = pow(g, (p - 1) // (2 * N), p)
    psi2 = pow(psi, 2, p)

    x = 1
    a_psi, b_psi = [], []
    for i in range(max(len(a), len(b))):
        a_psi.append(x * a[i] % p)
        b_psi.append(x * b[i] % p)
        x = x * psi % p

    C = [ai * bi % p for ai, bi in zip(ntt(a_psi, psi2), ntt(b_psi, psi2))]
    return intt(C, psi)


if __name__ == "__main__":
    trails = 100
    mul_times, schoolbook_times = [], []
    for _ in range(trails):
        a = [random.randrange(p) for _ in range(N)]
        b = [random.randrange(p) for _ in range(N)]

        start_time = time.perf_counter()
        c1 = negacyclic_mul(a, b)
        end_time = time.perf_counter()
        mul_times.append(end_time - start_time)

        start_time = time.perf_counter()
        c2 = schoolbook_negacyclic(a, b)
        end_time = time.perf_counter()
        schoolbook_times.append(end_time - start_time)

        assert c1 == c2
    print("[+] 負循環 NTT 乘法驗證通過")
    print(f"O(NlogN) 平均 {1000. * sum(mul_times) / trails:.3f} 毫秒")
    print(f"O(N^2) 平均 {1000. * sum(schoolbook_times) / trails:.3f} 毫秒")
