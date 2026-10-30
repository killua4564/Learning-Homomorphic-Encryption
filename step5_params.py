"""
步驟 5-2：改變參數，量化「深度、精度、速度、大小」的取捨

目的：
  步驟 4 在玩具上學到：想做更多次乘法 → 中間質數要更多 → q 更大 → 為了安全 N 要更大
  → 更慢、更大；或者縮小 Δ → 更不準。這一步用 TenSEAL 的真實參數把它量出來。

每組參數量測：
  - 建立 context + 產生 galois keys 的時間，以及公開金鑰（不含私鑰）序列化後的大小
  - 加密兩個向量、內積（.dot）、解密的時間，以及內積「平均每個 slot」的時間
  - 一個密文 serialize 之後的大小
  - 內積的相對誤差
  - 參數不合法時，記錄錯誤訊息

參數組（N、coeff_mod_bit_sizes、Δ 的 bit 數）：
  A  8192   [60, 40, 60]                 2^40   1 層（內積只需要 1 次乘法）
  B  8192   [60, 40, 40, 60]             2^40   2 層（步驟 5-1 的設定）
  C  8192   [60, 40, 40, 40, 60]         2^40   3 層，總和 240 > 218
  D  8192   [60, 30, 30, 60]             2^30   跟 B 一樣 2 層，只把 Δ 縮小（對照組）
  E  8192   [60, 30, 30, 30, 60]         2^30   縮小 Δ 來換到 3 層
  F  16384  [60, 40, 40, 60]             2^40   跟 B 一樣，只把 N 加倍
  G  16384  [60, 40 × 5, 60]             2^40   5 層
  H  16384  [60, 40 × 7, 60]             2^40   7 層（總和 400，接近上限 438）

結果（一次執行）：
  組  金鑰生成   公開金鑰   加密      內積       每 slot   密文      解密      相對誤差
  A   51 ms      18 MB      4.6 ms    5.1 ms     1.25 µs   230 KB    0.23 ms   1.4e-7
  B   89 ms      34 MB      4.4 ms    9.5 ms     2.32 µs   324 KB    0.46 ms   1.3e-7
  C   encryption parameters are not set correctly
  D   88 ms      30 MB      4.4 ms    10.7 ms    2.61 µs   275 KB    0.44 ms   4.4e-5
  E   144 ms     47 MB      4.9 ms    16.0 ms    3.90 µs   345 KB    0.62 ms   4.5e-5
  F   198 ms     74 MB      9.1 ms    21.6 ms    2.63 µs   643 KB    0.86 ms   1.4e-6
  G   630 ms     242 MB     13.9 ms   70.4 ms    8.59 µs   1.2 MB    2.3 ms    1.4e-6
  H   1051 ms    406 MB     48 ms     117 ms     14.3 µs   1.6 MB    3.8 ms    1.4e-6

觀察（只比較「只差一個條件」的兩組）：
  - Δ（B → D）：時間幾乎不變、密文略小，相對誤差大約 300 倍 → Δ 決定精度，不影響速度
    （成本看的是質數個數，不是質數多大：每個質數都用一個 64-bit 整數存、各做一次 NTT）
  - 層數（A → B、D → E、F → G → H）：每多一層所有成本都上升，誤差不變
    → 內積只用 1 次乘法，A 就是最佳參數，B 多留的一層讓時間和金鑰都將近翻倍
  - N（B → F）：單次內積 2.3 倍，但每 slot 只慢約 13%；相對誤差大 10 倍，
    因為 N = 16384 時質數必須 ≡ 1 (mod 32768)，離 2^40 較遠，scale 不一致的偏差較大
  - 安全上限（C）：總和超過 HE 標準表的 218 bit，SEAL 直接拒絕
  - 公開金鑰（大部分是 galois key）遠大於密文：A 為 18 MB（密文的 80 倍），H 達 406 MB；
    使用者第一次要上傳這份金鑰，之後每次查詢只需幾毫秒、幾百 KB

執行：python3 step5_params.py   （F 組會跑幾秒）
"""
import time

import numpy as np
import tenseal as ts

CONFIGS = [
    ("A", 8192, [60, 40, 60], 40),
    ("B", 8192, [60, 40, 40, 60], 40),
    ("C", 8192, [60, 40, 40, 40, 60], 40),
    ("D", 8192, [60, 30, 30, 60], 30),
    ("E", 8192, [60, 30, 30, 30, 60], 30),
    ("F", 16384, [60, 40, 40, 60], 40),
    ("G", 16384, [60, 40, 40, 40, 40, 40, 60], 40),
    ("H", 16384, [60, 40, 40, 40, 40, 40, 40, 40, 60], 40),
]


SPACES = (7, 7, 38, 8, 7, 17, 15, 15, 12, 17, 15, 15, 16)

def print_row(*args):
    spaces = SPACES
    if len(args) != len(SPACES):
        spaces = list(SPACES[:5]) + [sum(SPACES[5:]) + 7]

    assert len(args) == len(spaces)
    print("|" + "|".join(f"{arg:^{space}}" for arg, space in zip(args, spaces)) + "|")


if __name__ == "__main__":
    print_row("Group", "N", "coeff_mod_bit_sizes", "Δ bits", "depth", "generate key ms", "public key KB", "encryption ms", "dot ms", "average slot µs", "ciphertext KB", "decryption ms", "relative error")
    print_row(*["-" * i for i in SPACES])
    for name, n, coeff_mod_bit_sizes, scale_bits in CONFIGS:
        depth = len(coeff_mod_bit_sizes) - 2
        try:
            start_time_ns = time.perf_counter_ns()
            ctx = ts.context(ts.SCHEME_TYPE.CKKS, poly_modulus_degree=n, coeff_mod_bit_sizes=coeff_mod_bit_sizes)
            ctx.global_scale = 2 ** scale_bits
            ctx.generate_galois_keys()
            generate_key_ms = f"{(time.perf_counter_ns() - start_time_ns) / 1_000_000.:.3f} ms"
            public_key_kb = f"{len(ctx.serialize(save_secret_key=False)) / 1_024.:.2f} KB"

            rng = np.random.default_rng()
            a = rng.random(n // 2)
            b = rng.random(n // 2)
            expected = float(a @ b)

            start_time_ns = time.perf_counter_ns()
            enc_a = ts.ckks_vector(ctx, a)
            enc_b = ts.ckks_vector(ctx, b)
            encryption_ms = f"{(time.perf_counter_ns() - start_time_ns) / 1_000_000.:.3f} ms"
            ciphertext_kb = f"{len(enc_a.serialize()) / 1_024.:.2f} KB"

            start_time_ns = time.perf_counter_ns()
            enc_ab = enc_a.dot(enc_b)
            dot_ns = time.perf_counter_ns() - start_time_ns
            dot_ms = f"{dot_ns / 1_000_000.:.3f} ms"
            average_slot_us = f"{dot_ns / (n // 2) / 1_000.:.3f} µs"

            start_time_ns = time.perf_counter_ns()
            dec_ab = enc_ab.decrypt()[0]
            decryption_ms = f"{(time.perf_counter_ns() - start_time_ns) / 1_000_000.:.3f} ms"
            error = f"{100. * abs(expected - dec_ab) / expected:.2e}%"

            print_row(name, n, str(coeff_mod_bit_sizes), scale_bits, depth, generate_key_ms, public_key_kb, encryption_ms, dot_ms, average_slot_us, ciphertext_kb, decryption_ms, error)
        except Exception as error:
            print_row(name, n, str(coeff_mod_bit_sizes), scale_bits, depth, str(error))
