# 從 RSA 到格密碼：同態加密的玩具實作與雜訊實驗

從密碼學 CTF 熟悉的 RSA 出發，一步步實作部分同態加密（Paillier）、LWE 加密與格攻擊，再在玩具 LWE 上量測同態運算的雜訊成長，接著親手做出 CKKS 的 scale 與 rescale、再在多項式環上手刻一個 mini-CKKS（編碼、相乘、relinearization、旋轉），最後用 TenSEAL 量測真實參數下的速度、大小與精度。目的是用程式和實驗理解基於格的全同態加密（BFV／CKKS）背後的核心取捨：**雜訊預算、乘法深度、精度與參數安全性**。

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
| 4-1 | CKKS 的 scale | [`step4_scale.py`](step4_scale.py) | 實數乘上 Δ 再加密；誤差 ≈ 雜訊 / Δ；Δ·\|x\| 超過 q/2 會繞回 |
| 4-2 | 不 rescale 的乘法 | [`step4_mul.py`](step4_mul.py) | 倍數 Δ^k 一路長大；精度與深度在同一個 q 下互相拉扯 |
| 4-3 | rescale 與模數鏈 | [`step4_rescale.py`](step4_rescale.py) | 用質數鏈做 rescale、對齊 level 與 scale 計算 x² + x；單一玩具密文的 LLL 攻擊 |
| 5-0a | CKKS 編碼 | [`step5_0a_encode.py`](step5_0a_encode.py) | 把 N/2 個數字編碼成一個 X^N + 1 上的整數係數多項式；多項式加乘 = slot 逐格加乘 |
| 5-0b | RLWE 加解密 | [`step5_0b_rlwe.py`](step5_0b_rlwe.py) | 一個密文裝 4 個數字；sample extraction 只靠重排係數取得 slot 總和 |
| 5-0c | 多項式上的相乘與 rescale | [`step5_0c_mul.py`](step5_0c_mul.py) | 相乘後密文只多一個多項式；單一 slot 溢位會毀掉整個密文 |
| 5-0d | relinearization | [`step5_0d_relin.py`](step5_0d_relin.py) | 用「s² 的加密」把 3 個多項式壓回 2 個；特殊質數 P 的作用 |
| 5-0e | 旋轉與 galois key | [`step5_0e_rotate.py`](step5_0e_rotate.py) | X → X^g 讓 slot 換位置；用 galois key 換回鑰匙；做出 sum 與內積 |
| 5-1 | TenSEAL 內積 | [`step5_dot.py`](step5_dot.py) | 真實參數下的時間、大小、誤差；誤差主要來自 scale 不一致 |
| 5-2 | 參數取捨 | [`step5_params.py`](step5_params.py) | 改變 N、層數、Δ，量測金鑰、加解密、內積、密文與公開金鑰大小、誤差 |
| 5-3 | 做太深會怎樣 | [`step5_depth.py`](step5_depth.py) | 層數用完會報錯；數值溢位不會報錯、安靜地回傳錯誤結果 |
| 5-4 | SEAL 版 mini-CKKS | [`step5_seal.py`](step5_seal.py) | 把 mini-CKKS 改成 SEAL 的結構：NTT、RNS、每個質數拆開的 key switching；特殊質數 P 從 2⁹⁰ 降到 2²⁶ |

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

### 6. CKKS：scale、rescale 與參數取捨
在步驟 3 的玩具 LWE 上，把實數 x 乘上 scale Δ、四捨五入後放進相位（相位 = Δ·x + e），一步步做出 CKKS 的核心機制。

**scale 決定精度，q 決定上限**（q ≈ 2⁴⁰，x = π，每組 1000 次）：

| Δ | 2⁴ | 2⁸ | 2¹⁶ | 2²⁴ | 2³⁶ | 2³⁷ | 2³⁸ | 2⁴⁰ |
|---|---|---|---|---|---|---|---|---|
| 平均誤差 | 1.1×10⁻¹ | 6.8×10⁻³ | 2.6×10⁻⁵ | 1.1×10⁻⁷ | 2.6×10⁻¹¹ | 1.3×10⁻¹¹ | **4.00** | **3.00** |

- 誤差 ≈ 雜訊 / Δ：Δ 每大 2⁸ 倍，誤差就小 2⁸ 倍。雜訊大小由安全性決定，使用者只能透過 Δ 選擇精度。
- 必須 Δ·|x| < q/2。超過之後相位會繞回去，解出 x − k·(q/Δ)，變成完全錯誤的數（2³⁸、2⁴⁰ 兩欄）。CKKS 對雜訊是「精度逐漸變差」，對溢位則是「懸崖式出錯」。

**不 rescale 時，倍數會一路變成 Δ^(k+1)**（q ≈ 2⁶⁰，連續乘上 1.5）：

| Δ | 每次乘法吃掉約 | 預測可乘次數 | 實測第幾次壞掉 | 剛加密時的相對誤差 |
|---|---|---|---|---|
| 2⁵ | 5.6 bit | 9 | 10 | 約 2～10% |
| 2¹⁰ | 10.6 bit | 4 | 5 | 約 10⁻³ |
| 2²⁰ | 20.6 bit | 1 | 2 | 約 10⁻⁶ |

條件是 log₂(scale) + log₂|值| < log₂(q/2)。在同一個 q 下，**Δ 越小能乘越多次、但越不準**。

**rescale 與模數鏈**（q = 一個約 2²⁵ 的 base × 三個約 2²⁰ 的質數，Δ = 2²⁰）：
- rescale = 整個密文除以鏈尾的質數 p 並四捨五入：倍數 Δ² → Δ²/p ≈ Δ，模數 q → q/p。連續乘上新加密的密文時，能乘的次數等於中間質數的個數（3 次），跟總 bit 數相同、不 rescale 的版本一樣多。但如果是中間結果互相相乘（例如連續平方），不 rescale 時倍數的 bit 數每次翻倍（Δ → Δ² → Δ⁴），rescale 後每次只用掉 log Δ 個 bit，**能做的次數就多很多**（同樣的玩具參數連續平方：不 rescale 2 次、rescale 3 次；5-3 實驗 C 的參數：不 rescale 3 次、rescale 7 次）。
- rescale 的另一個用途是**對齊倍數**。計算 x² + x（正確值 3.75）時：不 rescale 直接相加，x 的倍數比 x² 小一個 Δ，解出來 x 整個消失（誤差 1.50）；rescale 後再相加，誤差降到 2.9×10⁻³。
- 剩下的誤差來自「質數 ≠ Δ」造成的 scale 不一致：質數離 2²⁰ 約 0.2% 時誤差 2.9×10⁻³，改用離 2²⁰ 只差 3～17 的質數後降到約 2×10⁻⁵，被加密雜訊蓋過。SEAL 也因此挑選非常接近 2^bit 數的質數，並在 scale 不一致時拒絕相加（`scale mismatch`）。

**對應到 TenSEAL 參數**：

| TenSEAL | 意義 |
|---|---|
| `global_scale` | Δ，決定精度 |
| `coeff_mod_bit_sizes` 第一個 | base：要裝得下最後的 Δ·\|x\| |
| 中間的數字 | 每次乘法 rescale 掉一個，個數 = 乘法深度 |
| 最後一個 | 特殊質數，只在 relinearization／旋轉時暫時使用，用來縮小換鑰匙帶來的雜訊 |
| 總和 | 必須 ≤ HE 標準 Table 1 的上限（三元私鑰、128-bit：N = 8192 → 218、16384 → 438、32768 → 881） |
| `auto_rescale`／`auto_mod_switch`／`auto_relin` | 本專案手動做的 rescale、level 對齊、鑰匙對齊，TenSEAL 預設自動處理 |

例如 Δ = 2⁴⁰ 想做 3 次乘法需要 `[60, 40, 40, 40, 60]` = 240 bit，超過 N = 8192 的 218：只能加大 N（變慢、密文變大，但一個密文能裝的數字也變多），或縮小 Δ（變不準）。**深度、精度、速度三者互相拉扯。**

**玩具參數不安全**：只拿到一個玩具密文 (c₀, c₁, c₂)（n = 2），用 Kannan embedding 建 4 維格，目標向量 (Δx + e, s₁W, s₂W, W) 約 2²²，而 Gaussian heuristic 約 2³⁶，LLL 每次都能找回私鑰。真實參數中一個 RLWE 密文相當於 N = 8192 組 LWE 樣本，格的維度高達數千，δ^d 太模糊而無法攻擊。

### 7. mini-CKKS：在多項式環上手刻 CKKS 的每個元件
用 N = 8（4 個 slot）的玩具參數，依序手刻 CKKS 的每個元件，再對應到 TenSEAL：

| 手刻的元件 | 做法 | TenSEAL／SEAL 的對應 |
|---|---|---|
| 編碼（5-0a） | 4 個 slot 值 + 共軛 → 解方程得到 8 個實數係數 → 乘 Δ 取整數 | `ts.ckks_vector(ctx, [...])` 的編碼 |
| 加解密（5-0b） | c = (−a·s + m + e, a)，相位 c₀ + c₁·s，乘法為 mod X^N + 1 的多項式乘法 | 加密、`.decrypt()` |
| 相乘（5-0c） | 把密文當成 s 的多項式相乘，part 數 = len₁ + len₂ − 1 | `enc1 * enc2` |
| rescale（5-0c） | 每個係數除以模數鏈尾端的質數並四捨五入 | `auto_rescale`、`coeff_mod_bit_sizes` 中間的質數 |
| relinearization（5-0d） | rlk = 「s² 的加密（倍數 P）」，用 c₂ 乘上它再除以 P | relin keys、`auto_relin`、最後一個特殊質數 |
| 旋轉（5-0e） | X → X^g 讓 slot 往左挪；galois key =「s(X^g) 的加密」把鑰匙換回 s | `generate_galois_keys()` |
| 加總／內積（5-0e） | 挪 N/4、…、1 格各一次再相加，log₂(N/2) 次旋轉 | `.sum()`、`.dot()` |

過程中的觀察：
- **只有乘法需要 mod X^N + 1**：slot 點都是 X^N + 1 的根，在這些點上 X^N 和 −1 的值相同。
- **slot 總和 = (N/2) × 訊息多項式的常數項**。因此只要重新排列密文係數（負循環矩陣的第一列），就能得到總和的 LWE 密文（sample extraction），完全不需要 galois key；代價是結果離開了「一個密文裝很多數字」的形式。
- **一個 slot 溢位，所有 slot 都壞**：溢位發生在多項式係數上，而每個係數都混著所有 slot。
- **金鑰裡的雜訊 e 是安全性、不是正確性的需要**：relinearization 金鑰去掉 e 仍能正確運算，但 b ≡ −a·s (mod P) 會變成無雜訊的線性方程，直接洩漏 s。
- **旋轉用的 g 不一定要是 5**：只要 g 的次方剛好走過 N/2 個 slot 點、且每對共軛只選到一個即可（N = 8 時 3 也可以）。CKKS 原始論文用 5，SEAL 用 3（[galois.h](https://github.com/microsoft/SEAL/blob/main/native/src/seal/util/galois.h) 的 `generator_ = 3`）。

**SEAL 版（5-4）**：5-0 把模數 q 當成一個大數，SEAL 不是這樣做的。5-4 把同一套 mini-CKKS 改成 SEAL 的結構，實驗結果與 5-0e 相同（加解密、連乘 3 次、旋轉、sum = 3.25、內積 = 0.75）：

| 元件 | 5-0 版 | 5-4 版（SEAL 的作法） |
|---|---|---|
| 係數 | 一個 mod q 的大數 | RNS：每個質數存一份餘數，只有解密時用 CRT 組回 |
| 多項式乘法 | 直接展開，O(N²) | 每個質數各做一次負循環 NTT，O(N log N) |
| rescale | 組回大數、除以 p、四捨五入 | c′ⱼ = (cⱼ − [c]_{p_l}) · p_l⁻¹ (mod pⱼ)，全程在小質數上完成 |
| key switching（relinearize、旋轉共用） | P ≈ 2⁹⁰（比整個 q 大），1 把金鑰 | 把 c₂ 拆成每個質數一個小整數，P ≈ 2²⁶（只比單一質數大），每個質數 1 把金鑰 |
| 編碼 | 解線性方程 | slot 依指數排成 ψ¹, ψ³, … 的順序，再用複數版 intt |
| 旋轉的介面 | 收 g 和一把 galois key | 跟 SEAL 的 `rotate_vector` 一樣：收「轉幾格」和整組 galois key |

過程中的觀察：
- **P 從 2⁹⁰ 降到 2²⁶，relinearize 前後誤差完全一樣**（1.14×10⁻⁵）。省下來的 bit 數可以拿去做乘法深度，因為 P 也算在 HE 標準表的安全上限裡。
- **RNS 的「每一份都小」不代表這個數小**：任何數的餘數都小於 p。如果每一份各乘一把金鑰，等於拿整個 c₂ 去乘，雜訊是 c₂·e/P。key switching 真正的分解是把每一份的餘數 [c₂]_{pᵢ} **拿出來當成獨立的小整數**，轉成所有質數（包括 P）的 RNS，再乘上第 i 把金鑰。雜訊就變成 Σ [c₂]_{pᵢ}·e/P。
- **代價是約 L² 的金鑰大小與計算量**：L 個資料質數需要 L 把金鑰，每把存成 L + 1 份。這是公開金鑰遠大於密文的原因之一（第 8 節：18 MB 對 230 KB）。另一個選擇是不拆、讓 P 比整個 q 大（5-0 版），但這樣 P 會吃掉安全上限。介於兩者之間的作法叫 hybrid key switching（[Han–Ki 2019](https://eprint.iacr.org/2019/688)）。
- **小係數的值要先用整數產生，再轉成 RNS**：s、e、m 都是如此。如果每一份各抽一個小數字，用 CRT 組回來會是一個大數。只有均勻亂數 a 可以在每個質數上各抽一份。
- **伺服器只拿得到密文和公開的轉換金鑰**：galois key 必須由使用者事先產生，缺了某個方向就轉不了（例如只有往左的金鑰時，就不能往右轉）。SEAL 預設產生往左、往右各 1、2、4、…、N/4 格加上共軛，共 2·log₂N − 1 把，N = 8192 時是 25 把。
- **仍與 SEAL 不同的地方**：SEAL 的密文一直存成 NTT 之後的值，乘法只要逐點相乘，galois 也只是重新排列；這裡每次乘法都做 ntt → 相乘 → intt。

### 8. TenSEAL 實測：參數取捨與失敗方式
**內積**（N = 8192、`[60, 40, 40, 60]`、Δ = 2⁴⁰）：加密兩個 4096 維向量約 4.5 ms、密文內積約 9.5 ms（numpy 明文約 1.2 µs，**約慢 8000 倍**）、一個密文約 323 KB、相對誤差約 1.3×10⁻⁷。這個誤差主要不是隨機雜訊，而是 rescale 除以的質數 p = 1099511480321 並非剛好 2⁴⁰：2⁴⁰/p − 1 = 1.341×10⁻⁷，與實測完全吻合，4096 個 slot 的偏差同向累加。

**參數取捨**（每組一次執行）：

| 組 | N | coeff_mod_bit_sizes | Δ | 層數 | 金鑰生成 | 公開金鑰 | 加密 | 內積 | 每 slot | 密文 | 解密 | 相對誤差 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| A | 8192 | [60, 40, 60] | 2⁴⁰ | 1 | 51 ms | 18 MB | 4.6 ms | 5.1 ms | 1.25 µs | 230 KB | 0.23 ms | 1.4×10⁻⁷ |
| B | 8192 | [60, 40, 40, 60] | 2⁴⁰ | 2 | 89 ms | 34 MB | 4.4 ms | 9.5 ms | 2.32 µs | 324 KB | 0.46 ms | 1.3×10⁻⁷ |
| C | 8192 | [60, 40, 40, 40, 60] | 2⁴⁰ | 3 | 參數被拒絕（總和 240 > 218） | | | | | | | |
| D | 8192 | [60, 30, 30, 60] | 2³⁰ | 2 | 88 ms | 30 MB | 4.4 ms | 10.7 ms | 2.61 µs | 275 KB | 0.44 ms | 4.4×10⁻⁵ |
| E | 8192 | [60, 30, 30, 30, 60] | 2³⁰ | 3 | 144 ms | 47 MB | 4.9 ms | 16.0 ms | 3.90 µs | 345 KB | 0.62 ms | 4.5×10⁻⁵ |
| F | 16384 | [60, 40, 40, 60] | 2⁴⁰ | 2 | 198 ms | 74 MB | 9.1 ms | 21.6 ms | 2.63 µs | 643 KB | 0.86 ms | 1.4×10⁻⁶ |
| G | 16384 | [60, 40×5, 60] | 2⁴⁰ | 5 | 630 ms | 242 MB | 13.9 ms | 70.4 ms | 8.59 µs | 1.2 MB | 2.3 ms | 1.4×10⁻⁶ |
| H | 16384 | [60, 40×7, 60] | 2⁴⁰ | 7 | 1051 ms | 406 MB | 48 ms | 117 ms | 14.3 µs | 1.6 MB | 3.8 ms | 1.4×10⁻⁶ |

只比較「只差一個條件」的兩組：
- **Δ 決定精度，不影響速度**（B → D）：時間幾乎相同，相對誤差約大 300 倍。成本取決於質數**個數**而非大小，因為每個質數都各占一個 64-bit 整數、各做一次 NTT。
- **層數決定成本，不影響精度**（A → B、D → E、F → G → H）：每多一層，時間、金鑰、密文都上升，誤差不變。內積只需 1 次乘法，因此 **A 才是內積的最佳參數**，B 多留的一層讓時間與金鑰都將近翻倍。
- **N 決定安全上限與一次能處理的數量**（B → F）：單次內積慢 2.3 倍，但每個 slot 只慢約 13%。相對誤差反而大 10 倍，因為 N = 16384 的質數必須 ≡ 1 (mod 32768)，離 2⁴⁰ 較遠，scale 不一致的偏差較大。
- **公開金鑰（主要是 galois key）遠大於密文**：A 組 18 MB（約密文的 80 倍），H 組達 406 MB。在「使用者上傳加密資料、伺服器計算」的情境中，使用者第一次要上傳這份金鑰，之後每次查詢只需數毫秒與數百 KB，因此**同一把金鑰重複使用越多次越划算**。

**挑參數的順序**：先數清楚需要幾次乘法決定層數 → 依精度需求決定 Δ → 查 HE 標準表選出能裝下這些質數的最小 N。

**兩種失敗方式**（把密文連續平方）：

| 情況 | 參數 | 結果 | TenSEAL 是否報錯 |
|---|---|---|---|
| 層數用完 | 1、2、7 層，x = 1.01 | 分別在第 2、3、8 次乘法停止 | ✅ `ValueError: scale out of bounds` |
| 數值溢位 | 7 層，x = 30 | 第 5 次起解出完全錯誤的數 | ❌ **沒有任何警告**，直到第 8 次層數用完才報錯 |

函式庫只記錄 scale 與 level，看不到明文，因此無法察覺數值溢位。溢位時機可以用 bit 數預測：Δ·值的 bit 數（40 + log₂30 × 2ᵏ）必須小於剩下的模數 − 1，k = 4 時 118.5 < 179，k = 5 時 197 > 139。**使用 CKKS 時，使用者必須事先掌握資料範圍**（例如先正規化輸入）。此外，相對誤差每平方一次約翻倍（(x(1+ε))² ≈ x²(1+2ε)），深度除了增加成本，也會累積誤差。

---

## 執行環境

- **Python 3.10 以上**
- **SageMath**（step2、step3_noise、step3_mult、step4 系列使用；開發時版本為 10.9）
- **pycryptodome**（step1）、**numpy**（step3_noise、step3_mult、step5 系列）
- **TenSEAL**（step5_dot、step5_params、step5_depth；開發時版本為 0.3.18）

```bash
python3 step1_paillier.py
sage -python step2_lwe.py
sage -python step2_lll_delta.py     # 約 40 秒
sage -python step3_noise.py         # 每個 k 跑 10000 次，需要數分鐘
sage -python step3_mult.py          # 約 30 秒
python3 step3_ntt.py
sage -python step4_scale.py
sage -python step4_mul.py
sage -python step4_rescale.py
python3 step5_0a_encode.py
python3 step5_0b_rlwe.py
python3 step5_0c_mul.py
python3 step5_0d_relin.py
python3 step5_0e_rotate.py
python3 step5_seal.py
python3 step5_dot.py
python3 step5_params.py
python3 step5_depth.py
```

`step3_noise.py` 裡的 `ZERO_COV` 開關，用來切換加密時隨機向量 x 的取法（{0, 1} 或 {−1, 0, 1}）。

---

## 限制

- 參數都是玩具規模，目的是讓雜訊與攻擊**看得見**，不代表真實方案的安全等級或效能。
- 步驟 3-2 的乘法沒有實作 relinearization 與 modulus switching，所以密文長度會指數成長，只能做十幾層。
- 步驟 4 的 CKKS 是單一數字的 LWE 版本，沒有多項式編碼（slots）、relinearization 與特殊質數，用來理解 scale 與 rescale，不是完整的 CKKS。
- 步驟 5-0 的 mini-CKKS 使用 N = 8、直接展開的多項式乘法與一個超過整個 q 的特殊質數，沒有 NTT 與 RNS（每個質數分開運算），結構上仍與 SEAL 不同。步驟 5-4 改成 SEAL 的結構，但仍是 N = 8，密文也沒有像 SEAL 一樣一直存成 NTT 形式。
- 步驟 5-2 每組只量一次，時間會有些微波動；公開金鑰大小以 `ctx.serialize(save_secret_key=False)` 計算。
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
- Z. Brakerski, C. Gentry, V. Vaikuntanathan, *Fully Homomorphic Encryption without Bootstrapping* (BGV). https://eprint.iacr.org/2011/277
- J. Fan, F. Vercauteren, *Somewhat Practical Fully Homomorphic Encryption* (BFV). https://eprint.iacr.org/2012/144
- J. H. Cheon, A. Kim, M. Kim, Y. Song, *Homomorphic Encryption for Arithmetic of Approximate Numbers* (CKKS), ASIACRYPT 2017. https://eprint.iacr.org/2016/421
- Microsoft SEAL. https://github.com/microsoft/SEAL
- OpenMined TenSEAL. https://github.com/OpenMined/TenSEAL
- Z. Brakerski, V. Vaikuntanathan, *Efficient Fully Homomorphic Encryption from (Standard) LWE*（relinearization）. https://eprint.iacr.org/2011/344
- C. Gentry, S. Halevi, N. P. Smart, *Homomorphic Evaluation of the AES Circuit*（特殊模數 P 的 key switching）. https://eprint.iacr.org/2012/099
- K. Han, D. Ki, *Better Bootstrapping for Approximate Homomorphic Encryption*（hybrid key switching）. https://eprint.iacr.org/2019/688
- V. Lyubashevsky, C. Peikert, O. Regev, *On Ideal Lattices and Learning with Errors over Rings*. https://eprint.iacr.org/2012/230
- Y. Elias, K. Lauter, E. Ozman, K. Stange, *Provably Weak Instances of Ring-LWE*, CRYPTO 2015. https://eprint.iacr.org/2015/106
- NIST FIPS 203 (ML-KEM). https://csrc.nist.gov/pubs/fips/203/final
- NIST FIPS 204 (ML-DSA). https://csrc.nist.gov/pubs/fips/204/final
- Homomorphic Encryption Standard v1.1. http://homomorphicencryption.org/wp-content/uploads/2018/11/HomomorphicEncryptionStandardv1.1.pdf
