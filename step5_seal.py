"""
步驟 5-4：SEAL 版 mini-CKKS —— RNS + NTT + 每個質數拆開的 key switching

目的：
  5-0 的 mini-CKKS 把模數 q = p0·p1·p2·p3 當成「一個大數」，SEAL 不是這樣做的。
  這一步把 5-0e 改造成跟 SEAL 一樣的結構，完成後每個函式都對得上 SEAL 的實作方式。

改造的 6 個部分（全部在這個檔案裡）：

  1. NTT：把步驟 3-3 的 ntt / intt / negacyclic_mul 改成 p、ψ 由參數傳入
     （PRIMES、PSIS 已經算好：每個質數都 ≡ 1 (mod 2N)，ψ 是它的 2N 次本原單位根）

  2. RNS：一個係數 x 不存 x 本身，而是存 (x mod p0, x mod p1, …)；一個多項式就變成
     「每個質數一份」的 list。加法、乘法在每個質數上各算各的（乘法用第 1 部分的 NTT），
     需要 x 本身時（只有解密時）再用 CRT 組回：
         Q = 所有質數的乘積，Q_i = Q / p_i
         x = Σ r_i · Q_i · (Q_i^(−1) mod p_i)   (mod Q)

  3. 加解密、相乘：跟 5-0e 一樣，只是多項式都換成 RNS 形式。
     密文仍是 (parts, scale, level)，每個 part 是「前 level 個質數的 RNS 多項式」。

  4. rescale：不用再組回大數。要除掉最後一個質數 p_l，對每個剩下的質數 p_j：
         c'_j = (c_j − [c]_{p_l}) · p_l^(−1)   (mod p_j)
     其中 [c]_{p_l} 是「mod p_l 那一份」中心化後的值。c − [c]_{p_l} 剛好能被 p_l 整除，
     所以這就是 4-3 的「除以 p 四捨五入」，只是全程在小質數上完成。丟掉 p_l 那一份，level 減 1。

  5. key switching（relinearize 和 rotate 共用）：
     5-0d 讓 P 比整個 q 還大，因為 c2 的係數跟 q 一樣大。SEAL 先把 c2 拆成每個質數的一份：
         c2 ≡ Σ_i [c2]_{p_i} · g_i   (mod Q)，g_i = Q_i · (Q_i^(−1) mod p_i)
     （g_i mod p_j 在 i = j 時是 1、否則是 0，這就是 CRT 的那組基底）
     每一份 [c2]_{p_i} 都小於 p_i，所以只要 P 比「單一質數」大就夠了。
     金鑰也跟著拆成每個 i 一把：ksk_i =「P · g_i · t 的加密」（模數 P·Q），t 是要換掉的鑰匙
         relinearize：t = s²；rotate：t = s(X^g)
     換鑰匙：在模數 P·q 底下算 Σ_i [c2]_{p_i} · ksk_i，再除以 P（用第 4 部分的公式除掉 P）
     → 得到 (d0, d1)，d0 + d1·s ≈ c2 · t。這個 P 也必須 ≡ 1 (mod 2N)，因為也要做 NTT。

  6. encode：改用 intt 的公式取代 solve_right（5-0a 第 4 題的推廣）

  另外，rotate 改成跟 SEAL 的 rotate_vector 一樣：收「往左轉幾格」和整組 galois key（gks）。

常數：
  PRIMES = [base ≈ 2^25, 三個 ≈ 2^20]，P ≈ 2^26（只比最大的質數大一點，不再需要 2^90）

結果（跟 5-0e 一樣的實驗）：
  - 加解密、密文相加：誤差約 5e-6 ～ 1e-5
  - z1·z2：mul（3 個 parts）1.14e-5 → relinearize（2 個 parts）1.14e-5 → rescale（level 4 → 3）1.20e-5
    → P 只有 2^26，relinearize 卻完全沒有增加誤差（5-0d 要 P ≈ 2^90 才做得到）
  - z3 連乘（mul → relinearize → rescale）：3 次都正確，第 4 次 level 只剩 base
  - 旋轉：往左 1、2 格正確；3 格 = 1 格 + 2 格；往右 1 格因為 gks 裡沒有 g = 13 的金鑰而失敗
  - sum_slots(z1) 每格都是 3.25；內積每格都是 0.75，誤差約 3e-5

觀察：
  - RNS 的「每一份都小」不代表這個數小（任何數的餘數都 < p）。如果每一份各乘一把金鑰，
    等於拿整個 c2 去乘，雜訊是 c2·e/P。key switching 要把每一份的餘數 [c2]_{p_i} 拿出來
    當成獨立的小整數，轉成所有質數（含 P）的 RNS，再乘第 i 把金鑰 → 雜訊變成 Σ [c2]_{p_i}·e/P
  - 代價：L 個資料質數 → L 把金鑰 × 每把 L + 1 份，金鑰大小和 key switching 的計算量都約 L²，
    這是 5-2 的公開金鑰（主要是 galois key）遠大於密文的原因之一。不拆的話 P 要比整個 q 大（5-0d），
    會吃掉安全上限的 bit 數；介於兩者之間的叫 hybrid key switching（https://eprint.iacr.org/2019/688）
  - 全程不組回大數，只有解密時用 CRT 組回（SEAL 的 decoder 也是）。rescale 和 key switching 的
    「除以質數」是同一個公式（rescale_part），[c]_{p_l} 必須中心化，否則變成無條件捨去
  - 小係數（s、e、m）要先用整數產生再轉 RNS；只有均勻亂數 a 可以在每個質數上各抽一份
    （每份各抽小數字，用 CRT 組回來會是大數）
  - 伺服器只拿得到密文、rlk、gks，拿不到 s，所以 galois key 要由使用者事先產生，缺了就轉不了。
    SEAL 預設產生往左、往右 1、2、4、…、N/4 格加上共軛，共 2·log2(N) − 1 把（N = 8192 時 25 把）
  - 仍與 SEAL 不同：SEAL 的密文一直存成 NTT 之後的值（乘法只要逐點相乘、galois 只是重新排列），
    這裡每次乘法都做 ntt → 相乘 → intt；SEAL 的 generator 是 3，這裡用 5

執行：python3 step5_seal.py
"""
import functools
import math
import random

N = 8
GENERATOR = 5
DELTA = 2**20
ZETA = math.e ** (1j * math.pi / N)
EXPONENTS = [pow(GENERATOR, j, 2 * N) for j in range(N // 2)]
SLOT_ROOTS = [ZETA ** exp for exp in EXPONENTS]

# ── 模數鏈 + 特殊質數（每個都 ≡ 1 (mod 2N)，可以做 NTT）────────────────────────
PRIMES = [33554273, 1048433, 1048273, 1048193]    # [base ≈ 2^25, 三個 ≈ 2^20]
PSIS = [20938261, 417443, 410851, 139356]         # 對應的 2N 次本原單位根 ψ
P = 67108753                                      # 特殊質數 ≈ 2^26，只要比單一質數大
P_PSI = 18112966
ALL_PRIMES = PRIMES + [P]
ALL_PSIS = PSIS + [P_PSI]


def encode(z, delta=DELTA):
    # 改用 intt 的公式，不再呼叫 solve_right
    target = [0] * N
    for exp, v in zip(EXPONENTS, z):
        target[(exp - 1) // 2] = v
        target[(2 * N - exp - 1) // 2] = v.conjugate()
    return [round(v.real * delta) for v in intt(target, ZETA)]


def decode(m, delta=DELTA):
    return [
        sum(root ** i * m[i] for i in range(N)).real / delta
        for root in SLOT_ROOTS
    ]


def ntt(a, w, p=None):
    n = len(a)
    if n == 1:
        return a

    if p is not None:
        w2 = pow(w, 2, p)
    else:
        w2 = w ** 2

    x = 1
    result = [0] * n
    evens, odds = ntt(a[::2], w2, p), ntt(a[1::2], w2, p)
    for i in range(n // 2):
        result[i] = evens[i] + x * odds[i]
        result[i + n // 2] = evens[i] - x * odds[i]
        x *= w

        if p is not None:
            result[i] %= p
            result[i + n // 2] %= p
            x %= p

    return result


def intt(c, psi, p=None):
    if p is not None:
        ipsi = pow(psi, -1, p)
        ipsi2 = pow(ipsi, 2, p)
        inv_n = pow(len(c), -1, p)
    else:
        ipsi = psi.conjugate()
        ipsi2 = ipsi ** 2
        inv_n = 1. / len(c)

    x = 1
    result = []
    for a in ntt(c, ipsi2, p):
        ci = inv_n * x * a
        x *= ipsi

        if p is not None:
            ci %= p
            x %= p

        result.append(ci)

    return result


def negacyclic_mul(a, b, p, psi):
    x = 1
    a_psi, b_psi = [], []
    for i in range(max(len(a), len(b))):
        a_psi.append(x * a[i] % p)
        b_psi.append(x * b[i] % p)
        x = x * psi % p

    psi2 = pow(psi, 2, p)
    c = [ai * bi % p for ai, bi in zip(ntt(a_psi, psi2, p), ntt(b_psi, psi2, p))]
    return intt(c, psi, p)


def crt(remainders, moduli):
    result = 0
    prod = functools.reduce(lambda a, b: a * b, moduli)
    for remainder, modulus in zip(remainders, moduli):
        m = prod // modulus
        result += remainder * pow(m, -1, modulus) * m
    return result % prod


def to_rns(poly, primes=PRIMES):
    return [[c % prime for c in poly] for prime in primes]


def from_rns(rns, primes=PRIMES):
    poly = []
    q = functools.reduce(lambda a, b: a * b, primes)
    for i in range(N):
        v = crt([rns_poly[i] for rns_poly in rns], primes)
        poly.append((v + q // 2) % q - q // 2)
    return poly


def rns_add(a, b, primes=PRIMES):
    return [
        [(ai + bi) % prime for ai, bi in zip(a_poly, b_poly)]
        for a_poly, b_poly, prime in zip(a, b, primes)
    ]


def rns_mul(a, b, primes=PRIMES, psis=PSIS):
    return [
        negacyclic_mul(a_poly, b_poly, prime, psi)
        for a_poly, b_poly, prime, psi in zip(a, b, primes, psis)
    ]


def keygen():
    return to_rns([random.randint(-1, 1) for _ in range(N)], ALL_PRIMES)


def encrypt(s, z):
    # a 直接在每個質數上各抽一份均勻亂數，m + e 用 to_rns 轉
    m = encode(z, DELTA)
    a = [[random.randint(0, prime - 1) for _ in range(N)] for prime in PRIMES]
    e = [random.randint(-3, 3) for _ in range(N)]

    neg_a = [[-v for v in a_poly] for a_poly in a]
    rns_m_e = to_rns([mi + ei for mi, ei in zip(m, e)])
    c = rns_add(rns_mul(neg_a, s), rns_m_e)
    return (c, a), DELTA, len(PRIMES)


def decrypt(s, ct):
    # 相位在 RNS 上算完，最後才 from_rns 組回、decode
    parts, scale, level = ct

    m, si = parts[0], s
    for part in parts[1:]:
        m = rns_add(m, rns_mul(part, si))
        si = rns_mul(si, s)

    return decode(from_rns(m, PRIMES[:level]), scale)


def mul(ct1, ct2):
    # level 不同時，較高那邊直接丟掉多出來的質數（RNS 的「降 level」就是少存幾份）
    (parts1, scale1, level1), (parts2, scale2, level2) = ct1, ct2
    level = min(level1, level2)

    parts = [[[0] * N for _ in range(level)] for _ in range(len(parts1) + len(parts2) - 1)]
    for i in range(len(parts1)):
        for j in range(len(parts2)):
            parts[i + j] = rns_add(parts[i + j], rns_mul(parts1[i], parts2[j]))

    return parts, scale1 * scale2, level


def rescale(ct):
    # c'_j = (c_j − [c]_{p_l}) · p_l^(−1) mod p_j，丟掉 p_l 那一份；scale / p_l，level − 1
    parts, scale, level = ct
    primes, p = PRIMES[:level - 1], PRIMES[level - 1]
    new_parts = [rescale_part(part, primes, p) for part in parts]
    return new_parts, scale / PRIMES[level - 1], level - 1


def rescale_part(part, primes, p):
    n = len(primes)
    neg_level = [(-v + p // 2) % p - p // 2 for v in part[-1]]
    rns_part = rns_add(part[:-1], [neg_level] * n, primes)
    for i in range(n):
        inv_p = pow(p, -1, primes[i])
        rns_part[i] = [c * inv_p % primes[i] for c in rns_part[i]]
    return rns_part


def ks_keygen(s, t):
    # s、t 都是 ALL_PRIMES 上的 RNS 多項式，t 是要換掉的鑰匙（relinearize：s²；rotate：s(X^g)）
    # 兩層迴圈，跑的 list 不一樣：
    #   外層「第幾把金鑰」：對每個 p_i ∈ PRIMES（不含 P）產生一把 ksk_i = (b_i, a_i)，共 len(PRIMES) 把
    #       Q = prod(PRIMES)（不含 P），Q_i = Q // p_i，g_i = Q_i · (Q_i^(−1) mod p_i)
    #       b_i = −a_i·s + e_i + P·g_i·t   (mod P·Q)
    #   內層「RNS 存幾份」：每把金鑰的 b_i、a_i 都是 ALL_PRIMES 上的 RNS 多項式
    #       （a_i 在每個質數上各抽一份均勻亂數；e_i 先抽整數再 to_rns；
    #         P·g_i 是常數，所以 P·g_i·t 可以每一份各自乘，rns_add 時再 mod）
    result = []
    prod = functools.reduce(lambda a, b: a * b, PRIMES)
    for pi in PRIMES:
        qi = prod // pi
        gi = qi * pow(qi, -1, pi)
        pgt = [[P * gi * v for v in t_poly] for t_poly in t]

        a = [[random.randint(0, prime - 1) for _ in range(N)] for prime in ALL_PRIMES]
        e = to_rns([random.randint(-3, 3) for _ in range(N)], ALL_PRIMES)

        neg_a = [[-v for v in a_poly] for a_poly in a]
        b = rns_add(rns_mul(neg_a, s, ALL_PRIMES, ALL_PSIS), rns_add(e, pgt, ALL_PRIMES), ALL_PRIMES)
        result.append((b, a))
    return result


def key_switch(c, keys, level):
    # c 是乘著「錯的鑰匙」的那個 part（RNS，前 level 個質數）
    # 1. 對每個 i < level：取 [c]_{p_i}（中心化），轉成 PRIMES[:level] + [P] 的 RNS
    # 2. 跟 keys[i] 的 b、a（只取前 level 個質數和 P 那一份）相乘，累加
    # 3. 用第 4 部分的公式除掉 P → 回傳 (d0, d1)
    primes = PRIMES[:level] + [P]
    psis = PSIS[:level] + [P_PSI]

    d0 = [[0] * N for _ in range(level + 1)]
    d1 = [[0] * N for _ in range(level + 1)]
    for ci, (bi, ai), pi in zip(c, keys, PRIMES[:level]):
        rns_ci = to_rns([(v + pi // 2) % pi - pi // 2 for v in ci], primes)
        d0 = rns_add(d0, rns_mul(rns_ci, bi[:level] + [bi[-1]], primes, psis), primes)
        d1 = rns_add(d1, rns_mul(rns_ci, ai[:level] + [ai[-1]], primes, psis), primes)

    return (rescale_part(d0, PRIMES[:level], P), rescale_part(d1, PRIMES[:level], P))


def relin_keygen(s):
    return ks_keygen(s, rns_mul(s, s, ALL_PRIMES, ALL_PSIS))


def relinearize(ct, rlk):
    # (c0, c1, c2) → (c0 + d0, c1 + d1)
    parts, scale, level = ct
    assert len(parts) == 3
    c0, c1, c2 = parts
    d0, d1 = key_switch(c2, rlk, level)
    return [rns_add(c0, d0), rns_add(c1, d1)], scale, level


def apply_galois(rns, g, primes=ALL_PRIMES):
    # 回傳 poly(X^g) mod (X^N + 1)，rns 的每一份各自處理（第 i 份 mod primes[i]）
    # 係數 poly[j] 移到 (j·g mod 2N)；若 ≥ N，就放到 (j·g mod 2N) − N 並變號
    # 變號寫成 p − v（跟 SEAL 一樣，結果不會是負數；v = 0 時會得到 p，之後的運算都會 mod，不影響）
    # 迴圈照 rns 的份數跑，ALL_PRIMES 的前 4 個就是 PRIMES：
    #   密文（level 份）只會用到前 level 個質數，s（5 份）會用到包括 P 在內的全部
    n = len(rns)
    new_rns = [[0] * N for _ in range(n)]
    for i in range(n):
        for j in range(N):
            k = j * g % (2 * N)
            if k < N:
                new_rns[i][k] = rns[i][j]
            else:
                new_rns[i][k - N] = primes[i] - rns[i][j]

    return new_rns


def galois_keygen(s, g):
    return ks_keygen(s, apply_galois(s, g))


def rotate(ct, steps, gks):
    # 跟 SEAL 的 rotate_vector 一樣：收「往左轉幾格」和整組 galois key，自己算 g、挑金鑰
    #   g = GENERATOR^steps mod 2N（steps 為負數就是往右轉，pow 會自動取反元素）
    #   gks 裡沒有這個 g 的金鑰就做不了（SEAL 也會丟出錯誤）→ 只能用有的金鑰組合，例如 3 格 = 1 格 + 2 格
    # c0、c1 都做 apply_galois 得到 (c0', c1')，這時 c1' 乘的是 s(X^g)（錯的鑰匙）
    # 整個 c1' 交給 key_switch 換成 (d0, d1)，所以新密文是 (c0' + d0, d1)，c1' 不留下來
    # scale 不會變，但要原封不動帶到新密文上
    parts, scale, level = ct
    assert len(parts) == 2
    g = pow(GENERATOR, steps, 2 * N)
    if g not in gks:
        raise KeyError(f"沒有往左轉 {steps} 格（g = {g}）的 galois key")

    c0, c1 = apply_galois(parts[0], g), apply_galois(parts[1], g)
    d0, d1 = key_switch(c1, gks[g], level)
    return [rns_add(c0, d0), d1], scale, level


def sum_slots(ct, gks):
    # 跟 5-0e 一樣：往左轉 1、2、…、N/4 格，每轉一次就加回自己，最後每一格都是總和
    # 伺服器只拿得到 ct 和 gks（不需要 s）；相加改用 rns_add
    _, scale, level = ct
    for k in range(int(math.log2(N)) - 1):
        rotated_ct = rotate(ct, 2 ** k, gks)
        ct =([rns_add(poly_i, poly_j) for poly_i, poly_j in zip(ct[0], rotated_ct[0])], scale, level)

    return ct


if __name__ == "__main__":
    s = keygen()
    rlk = relin_keygen(s)
    gks = {}
    for k in range(int(math.log2(N)) - 1):
        g = pow(GENERATOR, 2 ** k, 2 * N)
        gks[g] = galois_keygen(s, g)

    z1 = [1.5, -2.0, 3.25, 0.5]
    z2 = [2.0, 0.5, -1.0, 4.0]
    z3 = [1.1, -0.9, 1.2, 0.8]

    def show(name, expected, ct):
        # 印出：預期值、解密值、最大誤差，以及密文的 parts 個數、level
        decrypted = decrypt(s, ct)
        max_error = max(abs(x - y) for x, y in zip(expected, decrypted))
        print(f"{name}")
        print(f"  expected ={[round(v, 4) for v in expected]}")
        print(f"  decrypted={[round(v, 4) for v in decrypted]}")
        print(f"  len(parts)={len(ct[0])} level={ct[2]} 最大誤差 {max_error:.2e}")

    print(f"P ≈ 2^{P.bit_length()}（5-0d 是 2^90），gks 裡的 g = {sorted(gks)}")

    # 1. 加解密、密文相加
    print("\n1. 加解密、相加")
    enc1, enc2 = encrypt(s, z1), encrypt(s, z2)
    show("z1", z1, enc1)
    show("z2", z2, enc2)
    enc_add = ([rns_add(i, j) for i, j in zip(enc1[0], enc2[0])], DELTA, len(PRIMES))
    show("z1 + z2", [i + j for i, j in zip(z1, z2)], enc_add)

    # 2. z1 · z2：mul → relinearize → rescale，每一步都解密一次
    print("\n2. z1 · z2")
    z12 = [i * j for i, j in zip(z1, z2)]
    enc_mul = mul(enc1, enc2)
    show("mul（3 個 parts，用 1、s、s² 解密）", z12, enc_mul)
    enc_relin = relinearize(enc_mul, rlk)
    show("relinearize（P 只有 2^26，誤差應該不變）", z12, enc_relin)
    show("rescale（level 減 1）", z12, rescale(enc_relin))

    # 3. z3 連乘：每次乘上新加密的 z3 → relinearize → rescale，直到 level 用完
    print("\n3. z3 連乘")
    ct3, expected = encrypt(s, z3), z3
    for i in range(1, 5):
        if ct3[2] <= 1:
            print(f"第 {i} 次：level 只剩 base，沒有質數可以 rescale 了")
            break
        ct3 = rescale(relinearize(mul(ct3, encrypt(s, z3)), rlk))
        expected = [a * b for a, b in zip(expected, z3)]
        show(f"第 {i} 次", expected, ct3)

    # 4. 旋轉：gks 只有往左 1、2 格的金鑰；3 格要組合，-1 格（往右）沒有金鑰
    print("\n4. 旋轉")
    show("往左 1 格", z1[1:] + z1[:1], rotate(enc1, 1, gks))
    show("往左 2 格", z1[2:] + z1[:2], rotate(enc1, 2, gks))
    show("往左 3 格 = 1 格 + 2 格", z1[3:] + z1[:3], rotate(rotate(enc1, 1, gks), 2, gks))
    try:
        rotate(enc1, -1, gks)
    except KeyError as error:
        print(f"往右 1 格：{error}")

    # 5. sum 和內積（伺服器只拿得到密文、rlk、gks，拿不到 s）
    print("\n5. sum、內積")
    show("sum_slots(z1)：每一格都是總和", [sum(z1)] * (N // 2), sum_slots(enc1, gks))
    enc_dot = sum_slots(rescale(relinearize(mul(enc1, enc2), rlk)), gks)
    show("內積 z1·z2：mul → relinearize → rescale → sum_slots", [sum(z12)] * (N // 2), enc_dot)
