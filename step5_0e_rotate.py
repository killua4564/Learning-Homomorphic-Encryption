"""
步驟 5-0e：旋轉與 galois key —— 做出 sum 和內積

目的：
  一個密文裡的 4 個 slot 彼此獨立，相乘只會同位置對同位置。要把 4 個 slot 加起來，
  需要「旋轉」：把整串數字往左挪 k 格再跟自己相加（挪 2 格、再挪 1 格，第 0 格就是總和）。
  這一步要做出：明文上的旋轉 → 密文上的旋轉（會換掉鑰匙）→ 用 galois key 換回來 → sum → 內積。
  做完之後，TenSEAL 的 generate_galois_keys()、.sum()、.dot() 就都對得上了。

方案：
  【明文上的旋轉】把多項式裡的 X 全部換成 X^g（再 mod X^N + 1）：
      m(X) → m(X^g)
    slot j 原本是 m(ζ^(5^j))，換完之後變成 m(ζ^(5^j · g))。
    取 g = 5 時，5^j · 5 = 5^(j+1)，所以 slot j 拿到原本 slot j+1 的值 → 整串往左挪 1 格。
    取 g = 5^k mod 2N 就是往左挪 k 格。（slot 點和旋轉用的 g 要搭配；CKKS 原始論文用 5，SEAL 用 3。）
    實作：係數 m_i 原本乘的是 X^i，換完變成 X^(i·g)；i·g mod 2N 落在 [N, 2N) 時，
          因為 X^N = −1，要放到 (i·g − N) 的位置並變號。

  【密文上的旋轉】把 c0、c1 都換成 c0(X^g)、c1(X^g)：
      c0(X^g) + c1(X^g)·s(X^g) = 相位(X^g)
    → 解得開、slot 也確實轉了，但要用「換過的鑰匙」s(X^g) 才能解，不是原本的 s。
    （跟 5-0d 乘完之後要用 s² 一樣，都是「鑰匙對不上」。）

  【galois key：把 s(X^g) 換回 s】跟 5-0d 的 relinearization 一模一樣，只是 s² 換成 s(X^g)：
      gk_g = (b, a)，b = −a·s + e + P·s(X^g)   (mod P·Q)
      旋轉後的密文 (c0', c1') = (c0(X^g), c1(X^g))，c1' 乘的是 s(X^g)：
      用 c1' 去乘 gk_g 再除以 P → 得到 (d0, d1)，d0 + d1·s ≈ c1'·s(X^g)
      新密文 = (c0' + d0, d1)      ← 注意：c1' 整個被換掉，不是加上去
    每一個 g 都需要自己的一把 galois key。

結果：
  - 明文旋轉：g = 5 往左 1 格、g = 25 往左 2 格、g = 15（≡ −1 mod 16）每格變成共軛，實數不變
  - 只做 apply_galois 不換鑰匙：用 s 解密得到約 1e19 的亂數，用 s(X^g) 解密才正確
  - rotate（含 galois key）後直接用 s 解密：正確，且確實轉了
  - sum_slots：4 個 slot 都是 3.25；內積（mul → relinearize → rescale → sum_slots）4 個 slot 都是 0.75，
    誤差約 1e-5
  - 把 slot 點改成 3 的次方、旋轉用 g = 3，結果一樣能往左轉一格

觀察：
  - 為什麼要旋轉：slot 之間的運算只能同位置對同位置，旋轉是讓不同 slot 的值「碰到彼此」的唯一方法，
    內積的「加總」這一步需要它
  - 旋轉是對整個多項式做的，整串 slot 會一起移動，所以加總完每一格都是總和；
    用挪 N/4、N/8、…、1 格的方式，只需要 log2(N/2) 次旋轉（N = 8192 時 12 次，而不是 4095 次）
  - g 的選擇：只要 g 的次方剛好走過 N/2 個點、而且每對共軛只選到一個即可；CKKS 原始論文用 5，SEAL 用 3（galois.h 的 generator_ = 3）

執行：python3 step5_0e_rotate.py
"""
import copy
import itertools
import math
import random

# ── 5-0a 的編碼 ──────────────────────────────────────────────────────────────
N = 8
GENERATOR = 5
DELTA = 2**20
ZETA = math.e ** (1j * math.pi / N)
SLOT_ROOTS = [ZETA ** pow(GENERATOR, j, 2 * N) for j in range(N // 2)]
ALL_ROOTS = list(itertools.chain.from_iterable([[root, root.conjugate()] for root in SLOT_ROOTS]))


def encode(z, delta):
    roots = [[root ** i for i in range(N)] for root in ALL_ROOTS]
    target = list(itertools.chain.from_iterable([v, v.conjugate()] for v in z))
    return [round(v.real * delta) for v in solve_right(roots, target)]


def decode(m, delta):
    return [
        sum(root ** i * m[i] for i in range(N)).real / delta
        for root in SLOT_ROOTS
    ]

def solve_right(a, b):
    # return x which satisfy ax = b
    n = len(a)
    a, b = copy.deepcopy(a), copy.deepcopy(b)
    for i in range(n - 1):
        pivot = i
        for j in range(i + 1, n):
            if abs(a[pivot][i]) < abs(a[j][i]):
                pivot = j

        a[i], a[pivot] = a[pivot], a[i]
        b[i], b[pivot] = b[pivot], b[i]

        for j in range(i + 1, n):
            r = a[j][i] / a[i][i]
            for k in range(n):
                a[j][k] = a[j][k] - r * a[i][k]
            b[j] = b[j] - r * b[i]

    x = [0.] * n
    for i in range(n - 1, -1, -1):
        ax = sum(x[j] * a[i][j] for j in range(i + 1, n))
        x[i] = (b[i] - ax) / a[i][i]

    return x


# ── 模數鏈 + 特殊質數 ─────────────────────────────────────────────────────────
PRIMES = [33554393, 1048573, 1048571, 1048559]
P = 1237940039285380274899124191   # 特殊質數，約 2^90（比整個 Q ≈ 2^85 還大）


def modulus(level):
    q = 1
    for p in PRIMES[:level]:
        q *= p
    return q


# ── 多項式運算（5-0c）────────────────────────────────────────────────────────
def poly_add(a, b, q):
    return [(x + y) % q for x, y in itertools.zip_longest(a, b, fillvalue=0)]


def poly_mul(a, b, q):
    c = [0] * N
    for i in range(N):
        for j in range(N):
            if i + j < N:
                c[i + j] += a[i] * b[j]
            else:
                c[i + j - N] -= a[i] * b[j]
    return [v % q for v in c]


def poly_centered(a, q):
    return [(v + q // 2) % q - q // 2 for v in a]


# ── 加解密、相乘、rescale（5-0c）─────────────────────────────────────────────
def keygen():
    return [random.randint(-1, 1) for _ in range(N)]


def encrypt(s, z):
    level = len(PRIMES)
    q = modulus(level)
    m = encode(z, DELTA)
    a = [random.randint(0, q - 1) for _ in range(N)]
    e = [random.randint(-3, 3) for _ in range(N)]

    neg_a = [-v for v in a]
    c = poly_add(poly_mul(neg_a, s, q), poly_add(m, e, q), q)
    return [c, a], DELTA, level


def decrypt(s, ct):
    parts, scale, level = ct
    q = modulus(level)

    m, si = parts[0], s
    for part in parts[1:]:
        m = poly_add(m, poly_mul(part, si, q), q)
        si = poly_mul(si, s, q)

    return decode(poly_centered(m, q), scale)


def mul(ct1, ct2):
    (parts1, scale1, level1), (parts2, scale2, level2) = ct1, ct2
    level = min(level1, level2)
    q = modulus(level)

    parts = [[] for _ in range(len(parts1) + len(parts2) - 1)]
    for i in range(len(parts1)):
        for j in range(len(parts2)):
            parts[i + j] = poly_add(parts[i + j], poly_mul(parts1[i], parts2[j], q), q)

    return parts, scale1 * scale2, level


def rescale(ct):
    parts, scale, level = ct
    p = PRIMES[level - 1]
    q = modulus(level)

    new_parts = [
        [(v + p // 2) // p for v in poly_centered(part, q)]
        for part in parts
    ]
    return new_parts, scale / p, level - 1


# ── relinearization（5-0d）───────────────────────────────────────────────────
def relin_keygen(s):
    p, q = P, modulus(len(PRIMES))
    a = [random.randint(0, p * q - 1) for _ in range(N)]
    e = [random.randint(-3, 3) for _ in range(N)]

    neg_a = [-v for v in a]
    p_s2 = [p * v for v in poly_mul(s, s, q)]
    b = poly_add(poly_mul(neg_a, s, p * q), poly_add(e, p_s2, p * q), p * q)
    return (b, a)


def relinearize(ct, rlk):
    parts, scale, level = ct
    assert len(parts) == 3
    p, q = P, modulus(level)
    (c0, c1, c2), (b, a) = parts, rlk

    c2 = poly_centered(c2, q)
    def linear_c2(c):
        return [(v + p // 2) // p for v in poly_centered(poly_mul(c, c2, p * q), p * q)]

    c0, c1 = poly_add(c0, linear_c2(b), q), poly_add(c1, linear_c2(a), q)
    return [c0, c1], scale, level


# ── 這一步要寫的 ─────────────────────────────────────────────────────────────
def apply_galois(poly, g):
    # 回傳 poly(X^g) mod (X^N + 1)
    # 係數 poly[i] 移到 (i·g mod 2N)；若 ≥ N，就放到 (i·g mod 2N) − N 並變號
    # （不 mod q：s 這種小係數多項式也會用到；密文的係數之後 poly_add 會再 mod）
    new_poly = [0] * N
    for i in range(N):
        j = i * g % (2 * N)
        if j < N:
            new_poly[j] = poly[i]
        else:
            new_poly[j - N] = -poly[i]

    return new_poly


def galois_keygen(s, g):
    # 跟 relin_keygen 一樣，只是 P·s² 換成 P·s(X^g)
    p, q = P, modulus(len(PRIMES))
    a = [random.randint(0, p * q - 1) for _ in range(N)]
    e = [random.randint(-3, 3) for _ in range(N)]

    neg_a = [-v for v in a]
    p_sg = [p * v for v in apply_galois(s, g)]
    b = poly_add(poly_mul(neg_a, s, p * q), poly_add(e, p_sg, p * q), p * q)
    return (b, a)


def rotate(ct, g, gk):
    # ct 只能有 2 個 part
    # c0', c1' = apply_galois(c0, g), apply_galois(c1, g)
    # 用 c1'（中心化後）去乘 gk 的 b、a（mod P·q），中心化 → 除以 P 四捨五入，得到 d0、d1
    # 回傳 ([c0' + d0, d1], scale, level)
    parts, scale, level = ct
    assert len(parts) == 2
    p, q = P, modulus(level)
    (c0, c1), (b, a) = parts, gk

    c0, c1 = apply_galois(c0, g), apply_galois(c1, g)
    c1 = poly_centered(c1, q)
    def rotate_c1(c):
        return [(v + p // 2) // p for v in poly_centered(poly_mul(c, c1, p * q), p * q)]

    c0, c1 = poly_add(c0, rotate_c1(b), q), rotate_c1(a)
    return [c0, c1], scale, level


def sum_slots(s, ct):
    # gks = {g: galois key}
    # 往左挪 2 格（g = 25）加回自己，再往左挪 1 格（g = 5）加回自己
    # 提示：兩個密文相加 = 兩個 part 各自 poly_add（level、scale 要一樣）
    _, scale, level = ct
    q = modulus(level)

    for k in range(int(math.log2(N)) - 1):
        g = pow(GENERATOR, 2 ** k, 2 * N)
        gk = galois_keygen(s, g)
        rotated_ct = rotate(ct, g, gk)
        ct = ([poly_add(i, j, q) for i, j in zip(ct[0], rotated_ct[0])], scale, level)

    return ct


if __name__ == "__main__":
    s = keygen()
    rlk = relin_keygen(s)
    g_trails = (GENERATOR ** 0, GENERATOR ** 1, GENERATOR ** 2, 2 * N - 1)
    z1 = [1.5, -2.0, 3.25, 0.5]
    z2 = [2.0, 0.5, -1.0, 4.0]

    # 1. 明文上的旋轉：g = 5、25、15
    for g in g_trails:
        print(f"g={g}, decode={decode(apply_galois(encode(z1, DELTA), g), DELTA)}")

    # 2. 密文只做 apply_galois：分別用 s 和 apply_galois(s, g) 解密
    print()
    for g in g_trails:
        print(f"g={g}")
        enc_z1 = encrypt(s, z1)
        g_enc_z1 = ([apply_galois(poly, g) for poly in enc_z1[0]], *enc_z1[1:])
        print(f"用 s 解密 = {decrypt(s, g_enc_z1)}")
        print(f"用 galois(s) 解密 = {decrypt(apply_galois(s, g), g_enc_z1)}")

    # 3. rotate（含 galois key）後用 s 解密
    print()
    for g in g_trails:
        rotate_enc = rotate(encrypt(s, z1), g, galois_keygen(s, g))
        print(f"g={g} rotated decrypt={decrypt(s, rotate_enc)}")

    # 4. sum_slots 和內積
    print()
    dot = sum_slots(s, encrypt(s, z1))
    print(decrypt(s, dot))
    print(sum(z1))
    dot = sum_slots(s, rescale(relinearize(mul(encrypt(s, z1), encrypt(s, z2)), rlk)))
    print(decrypt(s, dot))
    print(sum(i * j for i, j in zip(z1, z2)))
