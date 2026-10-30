"""
步驟 5-1：用 TenSEAL（CKKS）加密兩個向量並做內積

目的：
  用真實參數跑一次完整流程，量測：誤差、加密時間、內積時間、密文大小。

跟步驟 4 玩具版的對應：
  ts.context(...)                    建立參數：N、模數鏈；跟你的 PRIMES 是同一件事
  ctx.global_scale = 2**40           Δ
  ts.ckks_vector(ctx, 一串數字)       加密：一個密文最多裝 N/2 個數字（slots），這裡 N=8192 → 4096 個
  enc1 * enc2                        逐個位置相乘；TenSEAL 會自動 relinearize + rescale
  enc1.dot(enc2)                     內積 = 逐個位置相乘，再把所有位置加總
  .decrypt()                         解密，回傳 list
  .serialize()                       轉成 bytes，len() 就是傳輸／儲存時的大小

為什麼內積需要 galois keys：
  一個密文裡的 4096 個 slot 彼此獨立，相乘只會「同位置對同位置」。
  要把 4096 個位置加起來，需要「旋轉」：把整串數字往左挪 k 格，再跟原本的相加。
  挪 2048、1024、…、1 格各做一次（共 12 次），所有位置就加到同一格了。
  旋轉跟 relinearization 一樣是「換鑰匙」的運算，需要事先產生的公開金鑰，
  也就是 ctx.generate_galois_keys() 產生的那組（也會用到最後那個特殊質數）。

結果（N = 8192、[60, 40, 40, 60]、Δ = 2^40，三次執行）：
  加密兩個向量    約 4～5 ms
  密文內積        約 9.5 ms（numpy 明文內積約 1.2 µs，約慢 8000 倍）
  內積誤差        約 1.35e-4，相對誤差約 1.3e-7
  一個密文        約 323 KB（估算 2 × 8192 × 3 個質數 × 8 byte ≈ 384 KB，序列化時有壓縮）

觀察：
  - 相對誤差 1.3e-7 主要不是隨機雜訊，而是 scale 不一致：rescale 除以的質數 p = 1099511480321
    不是剛好 2^40，2^40 / p − 1 = 1.341e-7，每個 slot 都往同一個方向偏，4096 個加總後一起放大
    （SEAL 挑質數的規則見 native/src/seal/util/numth.cpp 的 get_primes：從 2^bit 往下每次減 2N，
     找到質數就收）
  - 把 generate_galois_keys() 拿掉，dot 會報錯：內積的加總需要旋轉，旋轉需要 galois key

安全規則：只使用隨機產生的資料。

執行：python3 step5_dot.py
"""
import time

import numpy as np
import tenseal as ts

N = 8192
COEFF_MOD_BIT_SIZES = [60, 40, 40, 60]
SCALE_BITS = 40
SLOTS = N // 2


def make_context():
    # 1. ctx = ts.context(ts.SCHEME_TYPE.CKKS, poly_modulus_degree=N, coeff_mod_bit_sizes=COEFF_MOD_BIT_SIZES)
    # 2. 設定 ctx.global_scale
    # 3. ctx.generate_galois_keys()
    ctx = ts.context(ts.SCHEME_TYPE.CKKS, poly_modulus_degree=N, coeff_mod_bit_sizes=COEFF_MOD_BIT_SIZES)
    ctx.global_scale = 2**40
    ctx.generate_galois_keys()
    return ctx


if __name__ == "__main__":
    ctx = make_context()

    # 隨機資料：兩個長度 4096、介於 0～1 的向量
    rng = np.random.default_rng()
    a = rng.random(SLOTS)
    b = rng.random(SLOTS)
    expected = float(a @ b)          # 明文內積，用來對答案

    # 加密 a、b，用 time.perf_counter() 量加密時間（兩個一起算或各自算都可以）
    start_time = time.perf_counter()
    enc_a = ts.ckks_vector(ctx, a)
    enc_b = ts.ckks_vector(ctx, b)
    end_time = time.perf_counter()
    print(f"加密 a 和 b 用了 {1000. * (end_time - start_time):.3f} 毫秒")

    # 在密文上做內積，量運算時間
    start_time = time.perf_counter()
    enc_ab = enc_a.dot(enc_b)
    end_time = time.perf_counter()
    print(f"密文做內積用了 {1000. * (end_time - start_time):.3f} 毫秒")

    # 解密（結果是只有一個元素的 list），算出絕對誤差與相對誤差
    ab = enc_ab.decrypt()[0]
    print(f"內積誤差 {abs(expected - ab):.2e}")

    # 量一個密文 serialize() 之後的大小（KB）
    print(f"密文序列化的大小 {len(enc_a.serialize()) / 1024:.2f} KB")

    # 把上面的數字印出來
    print(f"{ab=} {expected=}")
