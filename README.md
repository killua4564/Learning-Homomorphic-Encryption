# 從 RSA 到格密碼：同態加密的玩具實作與雜訊實驗

從密碼學 CTF 熟悉的 RSA 出發，一步步實作部分同態加密（Paillier）、LWE 加密與格攻擊，再在玩具 LWE 上量測同態運算的雜訊成長。目的是用程式和實驗理解基於格的全同態加密（BFV／CKKS）背後的核心取捨：**雜訊預算、乘法深度與參數安全性**。

> ⚠️ 所有實作都是**教學用的玩具版本**，參數刻意取得很小，方便觀察與攻擊，**不具備實際安全性**，請勿用於真實資料。

---

## 學習路線

| 步驟 | 主題 | 檔案 | 重點 |
|---|---|---|---|
| 1 | 部分同態 | [`step1_paillier.py`](step1_paillier.py) | RSA 乘法同態 → Paillier 加法同態；為什麼只能做一種運算、為什麼量子電腦能破 |
| 2 | LWE 與格攻擊 | [`step2_lwe.py`](step2_lwe.py) | Regev LWE 位元加密；Kannan embedding + LLL 從公鑰回推私鑰 |
| 2 | LLL 的能力上限 | [`step2_lll_delta.py`](step2_lll_delta.py) | 實測 LLL 的 root Hermite factor δ |
| 3-1 | 加法的雜訊成長 | [`step3_noise.py`](step3_noise.py) | 密文連加的雜訊分布；共用金鑰造成的共變異數 |
| 3-2 | 乘法與雜訊預算 | [`step3_mult.py`](step3_mult.py) | 用 outer product 做密文乘法；乘法深度與 log q 的關係 |
| 3-3 | RLWE 與 NTT | [`step3_ntt.py`](step3_ntt.py) | X^N + 1 上的負循環 NTT，O(N log N) 多項式乘法 |

---

## 主要發現

### 1. 部分同態只能做一種運算
Paillier 把明文放在指數上（g = n + 1，利用 (1 + n)^m ≡ 1 + mn mod n²），所以「密文相乘 = 明文相加」。若要做到明文相乘，等於要從 g^a、g^b 算出 g^(ab)，這是 Diffie–Hellman 形狀的難題。RSA 和 Paillier 都建立在因數分解上，會被 Shor 演算法（量子求週期）破解。

### 2. LWE 的困難來自「雜訊 + mod q」
沒有雜訊時，b = A·s 用高斯消去就能解；有雜訊但沒有 mod 時，可以用最小平方法近似。兩者同時存在，才變成困難問題。

把 LWE 轉成格的短向量問題後，用 LLL 攻擊（每組 5 次）：

| m（樣本數） | 攻擊成功的 n | 開始失敗的 n |
|---|---|---|
| 40 | ≤ 26 | 28 |
| 80 | ≤ 34 | 40 |

能不能攻擊成功，取決於目標短向量 (e, ±1) 相對於格的 Gaussian heuristic 有多「突出」（間隙比），以及 LLL 的模糊程度 δ^d。實測 δ 在 d = 40～160 時約為 1.015～1.021，與文獻中的 1.0219 一致。

由此得到參數方向：**n、雜訊 B 變大 → 較難攻擊；q 變大 → 較易攻擊**。但解密正確需要 q/B 大，安全性需要 q/B 小，只能靠加大 n 換取空間，這就是同態加密參數取捨的核心。

### 3. 加法雜訊：共用金鑰讓雜訊長得比預期快
完整數據（含失敗率與平均偏移）可執行 `step3_noise.py` 重現。

把 k 個 LWE 密文相加（n = 26, q = 97, m = 40，每個 k 跑 10000 次）：

| k | 實測標準差 | 若各自獨立（√k） | 考慮共變異數 √((20/3)·k(k+1)) | x 改為以 0 為中心後 |
|---|---|---|---|---|
| 1 | 3.65 | 3.65 | 3.65 | 4.19（預測 4.22） |
| 2 | 6.37 | 5.16 | 6.32 | 6.00（預測 5.96） |
| 5 | 14.03 | 8.16 | 14.14 | 9.38（預測 9.43） |
| 10 | 23.97 | 11.55 | 27.08* | 13.46（預測 13.33） |

\* k ≥ 10 時雜訊開始被 mod q 折回，接近飽和（約 28，即 [−48, 48] 上均勻分布的標準差），此時失敗率約 50%，密文不再含有任何明文資訊。

- 同一把金鑰加密的密文共用同一組 e，彼此有**共變異數**，雜訊大約跟 **k** 成長，而不是 √k。
- 把加密用的隨機向量 x 從 {0, 1} 改為以 0 為中心的 {−1, 0, 1}，共變異數就消失，雜訊回到 √k 成長。k = 10 時失敗率從 36% 降到 7.5%。ML-KEM 加密時採用中心化的分布，也是同樣的考量。
- q = 97 是奇數，2·⌊q/2⌋ ≡ −1，因此每兩個 1 相加會造成 −1 的偏移，平均偏移為 −(k − 1)/4，與實測吻合。

### 4. 乘法雜訊：雜訊預算與乘法深度
用 BGV 式的低位編碼（相位 = m + t·e）搭配 outer product 做出不需私鑰的密文乘法，連續乘上新鮮密文：

| q | log₂ q | 預測可撐次數 (log₂ q − 1) / 1.52 | 實測平均第幾次出錯 |
|---|---|---|---|
| 2⁸ + 1 | 8 | 4.6 | 4.3 |
| 2¹⁴ + 1 | 14 | 8.6 | 8.3 |
| 2²⁰ + 1 | 20 | 12.5 | 10.6（上限 11 次，被低估） |

- 每乘一次，雜訊平均增加 E[log₂|1 + 2e|] ≈ 1.52 bit。log₂(q/2) 就是**雜訊預算**，能做的乘法次數（深度）與 log q 成正比。
- 密文長度每乘一次就變成 (n + 1) 倍，這就是真實方案需要 **relinearization** 的原因。
- **q 必須與明文模數 t 互質**：若 t | q，把密文 mod t 會讓雜訊 t·e 消失，LWE 退化成無雜訊的線性方程組，可直接解出私鑰。

### 5. RLWE 與 NTT
- Z[X]/(X^N + 1) 上的乘法等於「負循環矩陣 × 向量」：一組 RLWE 樣本相當於 N 組結構化的 LWE 樣本。金鑰從 n² 個數降到 n 個，一個密文也能裝 N 個數。
- 手寫遞迴負循環 NTT：N = 256 時約 0.55 毫秒，直接展開約 3.2 毫秒。
- X^N + 1 的根是 2N 次本原單位根 ψ 的奇數次方，因此 NTT 需要 **q ≡ 1 (mod 2N)**，這正是 SEAL／TenSEAL 對 `coeff_mod_bit_sizes` 中每個質數的要求。
- 不能改用 X^N − 1：X = 1 是它的根，代入後 b(1) = a(1)·s(1) + e(1)，雜訊仍然很小，問題退化成一維 LWE。

---

## 執行環境

- **Python 3.10 以上**
- **SageMath**（step2、step3_noise、step3_mult 使用；開發時版本為 10.9）
- **pycryptodome**（step1）、**numpy**（step3_noise、step3_mult）

```bash
python3 step1_paillier.py
sage -python step2_lwe.py
sage -python step2_lll_delta.py     # 約 40 秒
sage -python step3_noise.py         # 每個 k 跑 10000 次，需要數分鐘
sage -python step3_mult.py          # 約 30 秒
python3 step3_ntt.py
```

`step3_noise.py` 裡的 `ZERO_COV` 開關，用來切換加密時隨機向量 x 的取法（{0, 1} 或 {−1, 0, 1}）。

---

## 限制

- 參數都是玩具規模，目的是讓雜訊與攻擊**看得見**，不代表真實方案的安全等級或效能。
- 步驟 3-2 的乘法沒有實作 relinearization 與 modulus switching，所以密文長度會指數成長，只能做十幾層。
- 雜訊量測使用 `bit_length` 取整數，每步增加的 bit 數會比理論值略低。
- 實驗次數有限（LLL 攻擊每組 5 次、乘法每組 50 次），數字是趨勢，不是精確門檻。

---

## 參考資料

- P. Paillier, *Public-Key Cryptosystems Based on Composite Degree Residuosity Classes*, EUROCRYPT 1999. https://link.springer.com/chapter/10.1007/3-540-48910-X_16
- P. W. Shor, *Polynomial-Time Algorithms for Prime Factorization and Discrete Logarithms on a Quantum Computer*. https://arxiv.org/abs/quant-ph/9508027
- O. Regev, *On Lattices, Learning with Errors, Random Linear Codes, and Cryptography*. https://cims.nyu.edu/~regev/papers/qcrypto.pdf
- O. Regev, *The Learning with Errors Problem* (survey). https://cims.nyu.edu/~regev/papers/lwesurvey.pdf
- N. Gama, P. Q. Nguyen, *Predicting Lattice Reduction*, EUROCRYPT 2008. https://www.iacr.org/archive/eurocrypt2008/49650031/49650031.pdf
- M. Albrecht et al., *On the Concrete Hardness of Learning with Errors*. https://eprint.iacr.org/2015/046
- V. Lyubashevsky, C. Peikert, O. Regev, *On Ideal Lattices and Learning with Errors over Rings*. https://eprint.iacr.org/2012/230
- Y. Elias, K. Lauter, E. Ozman, K. Stange, *Provably Weak Instances of Ring-LWE*, CRYPTO 2015. https://eprint.iacr.org/2015/106
- NIST FIPS 203 (ML-KEM). https://csrc.nist.gov/pubs/fips/203/final
- NIST FIPS 204 (ML-DSA). https://csrc.nist.gov/pubs/fips/204/final
- Homomorphic Encryption Standard v1.1. http://homomorphicencryption.org/wp-content/uploads/2018/11/HomomorphicEncryptionStandardv1.1.pdf
