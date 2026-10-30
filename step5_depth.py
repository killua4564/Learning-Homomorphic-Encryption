"""
步驟 5-3：故意做太深 —— 在哪裡失敗？怎麼失敗？

目的：
  步驟 4、5-0 在玩具上看過兩種壞法：
    (1) 質數用完（level 用完）：沒有空間再 rescale
    (2) 數值太大：Δ·|x| 超過「目前剩下的模數」，相位繞回去，解出錯的數（4-1、5-0c 的 z1）
  這一步在 TenSEAL 上重現這兩種失敗，記錄「第幾次乘法失敗」和「錯誤訊息」。

實驗：
  每組參數把同一個密文連續平方（enc = enc * enc），每一步解密、算相對誤差，
  直到丟出例外為止。

  實驗 1：數值維持在 1 附近（x = 1.01），只看層數用完
    A  8192   [60, 40, 60]          1 層
    B  8192   [60, 40, 40, 60]      2 層
    C  16384  [60, 40 × 7, 60]      7 層
  實驗 2：數值會長很大（x = 30，平方 k 次後是 30^(2^k)）
    D  16384  [60, 40 × 7, 60]      7 層

結果：
  實驗 1：A 在第 2 次、B 在第 3 次、C 在第 8 次乘法丟出 ValueError: scale out of bounds
          → 能做的乘法次數 = 中間質數的個數（1、2、7）
  實驗 2：D 前 4 次正確，第 5 次起解出完全錯誤的數（相對誤差 100%），但沒有任何例外，
          一直到第 8 次層數用完才丟出 scale out of bounds

觀察：
  - 兩種失敗，TenSEAL 只抓得到一種：
      層數用完 → 會報錯（函式庫自己記錄 scale 和 level）
      數值溢位 → 不會報錯（函式庫看不到明文，不知道值有多大），安靜地回傳錯的數
    → 使用 CKKS 時，使用者必須事先掌握資料範圍（例如先把輸入正規化到小範圍）
  - 用 bit 數預測 D 的溢位：Δ·值 的 bit 數 = 40 + log2(30)·2^k，必須小於 剩下的模數 − 1
      k = 4：40 + 78.5 = 118.5 < 60 + 3×40 − 1 = 179 → 正確
      k = 5：40 + 157  = 197   > 60 + 2×40 − 1 = 139 → 溢位
  - 相對誤差每平方一次大約翻倍（(x(1 + ε))² ≈ x²(1 + 2ε)）：C 組從 1.4e-6 一路長到 3.2e-4，
    深度除了增加成本，也會累積誤差

執行：python3 step5_depth.py
"""
import math

import numpy as np
import tenseal as ts

EXPERIMENTS = [
    ("A", 8192, [60, 40, 60], 1.01),
    ("B", 8192, [60, 40, 40, 60], 1.01),
    ("C", 16384, [60, 40, 40, 40, 40, 40, 40, 40, 60], 1.01),
    ("D", 16384, [60, 40, 40, 40, 40, 40, 40, 40, 60], 30.0),
]
SCALE_BITS = 40
MAX_MULS = 10


if __name__ == "__main__":
    for name, n, coeff_mod_bit_sizes, x in EXPERIMENTS:
        # 1. 建立 context，global_scale = 2 ** SCALE_BITS
        # 2. 加密一個全部都是 x 的向量（長度 n // 2）
        # 3. 最多做 MAX_MULS 次：enc = enc * enc，expected 也平方
        #    每一步解密，印出：第幾次、正確值、解出的值、相對誤差、
        #    以及 log2(Δ·正確值) 和「還剩下的模數 bit 數」
        # 4. 丟出例外時，印出「第幾次乘法」和錯誤訊息，結束這一組

        print(f"{name=} {n=} {coeff_mod_bit_sizes=} {x=}")
        ctx = ts.context(ts.SCHEME_TYPE.CKKS, poly_modulus_degree=n, coeff_mod_bit_sizes=coeff_mod_bit_sizes)
        ctx.global_scale = 2 ** SCALE_BITS
        ctx.generate_galois_keys()

        try:
            expected = x
            enc = ts.ckks_vector(ctx, [x] * (n // 2))
            dec = enc.decrypt()[0]
            error_rate = 100. * abs(dec - expected) / expected
            print(f"第 {0} 次 expected={expected:.4f} decrypt={dec:.4f} error_rate={error_rate:.2e}%")
            for i in range(1, MAX_MULS + 1):
                expected = expected * expected
                enc = enc * enc
                dec = enc.decrypt()[0]
                error_rate = 100. * abs(dec - expected) / expected
                print(f"第 {i} 次 expected={expected:.4f} decrypt={dec:.4f} error_rate={error_rate:.2e}%")
        except Exception as error:
            print(error)
