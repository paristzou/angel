beverage_edgeR_practice.R
================
paris
2026-10-25

``` r
# ==============================================================================
# 連鎖手搖飲銷量大數據分析實作 (Beverage Store Sales Pipeline with edgeR)
# 
# 授課主題：以手搖飲門市銷量理解 RNA-seq 差異表達分析 (edgeR QL Workflow)
# 對照概念：
#   - 基因 (Gene)            <--> 飲料單品 (Beverage Item, 如波霸厚鮮奶、四季春)
#   - 定序計數 (Read Count)  <--> 當日銷售杯數 (Daily Cups Sold)
#   - 定序深度 (Library Size)<--> 門市當日總營收杯數 (Total Store Daily Volume)
#   - 實驗分組 (Group)       <--> 門市類型 x 活動 (Campus.Normal, Campus.Promo...)
#   - 基因集 (Gene Set)      <--> 飲料分類 (Category, 如季節冰沙、原葉純茶)
# ==============================================================================
```

``` r
# 步驟 0. 安裝與載入必備套件 (Setup Packages)
```

``` r
# 若尚未安裝 Bioconductor 或 edgeR，請取消下行註解執行：
# if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# BiocManager::install(c("edgeR", "limma"))
# install.packages(c("readxl", "gplots", "RColorBrewer"))

suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
  library(readxl)
  library(gplots)
  library(RColorBrewer)
})

cat(">>> 成功載入分析套件！開始執行手搖飲銷售數據分析 Pipeline...\n")
```

    ## >>> 成功載入分析套件！開始執行手搖飲銷售數據分析 Pipeline...

``` r
# 步驟 1. 讀取 Excel 資料並建立 DGEList 核心容器 (Data Import & DGEList)
```

``` r
excel_path <- "input.xlsx"

# 1.1 讀取樣本實驗設計表 (Store_Targets)
targets <- as.data.frame(read_excel(excel_path, sheet = "Store_Targets"))
rownames(targets) <- targets$SampleID
cat("\n--- [1.1 門市樣本實驗設計表] ---\n")
```

    ## 
    ## --- [1.1 門市樣本實驗設計表] ---

``` r
print(head(targets))
```

    ##                          SampleID StoreType Status         Group
    ## Campus_Normal_R1 Campus_Normal_R1    Campus Normal Campus.Normal
    ## Campus_Normal_R2 Campus_Normal_R2    Campus Normal Campus.Normal
    ## Campus_Normal_R3 Campus_Normal_R3    Campus Normal Campus.Normal
    ## Campus_Promo_R1   Campus_Promo_R1    Campus  Promo  Campus.Promo
    ## Campus_Promo_R2   Campus_Promo_R2    Campus  Promo  Campus.Promo
    ## Campus_Promo_R3   Campus_Promo_R3    Campus  Promo  Campus.Promo
    ##                  DayOfWeek Weather Replicate
    ## Campus_Normal_R1       Fri   Sunny         1
    ## Campus_Normal_R2       Fri   Sunny         2
    ## Campus_Normal_R3       Tue  Cloudy         3
    ## Campus_Promo_R1        Sat   Sunny         1
    ## Campus_Promo_R2        Sat   Sunny         2
    ## Campus_Promo_R3        Sun  Cloudy         3

``` r
# 1.2 讀取銷售杯數矩陣 (Sales_Counts)
counts_df <- as.data.frame(read_excel(excel_path, sheet = "Sales_Counts"))
rownames(counts_df) <- counts_df$ItemID
counts_mat <- as.matrix(counts_df[, -1]) # 移除 ItemID 欄位轉為純數值矩陣
storage.mode(counts_mat) <- "integer"

# 1.3 讀取飲品資訊表 (Beverage_Info, 相當於基因註解 Gene Annotation)
item_info <- as.data.frame(read_excel(excel_path, sheet = "Beverage_Info"))
rownames(item_info) <- item_info$ItemID

# 1.4 定義實驗分組因子 (Group Factor: 門市類型.促銷狀態)
group <- factor(targets$Group)
cat("\n--- [1.4 門市實驗分組 (Levels)] ---\n")
```

    ## 
    ## --- [1.4 門市實驗分組 (Levels)] ---

``` r
print(table(group))
```

    ## group
    ## Campus.Normal  Campus.Promo   Mall.Normal    Mall.Promo 
    ##             3             3             3             3

``` r
# 1.5 封裝為 edgeR 核心容器 DGEList
y <- DGEList(counts = counts_mat, group = group, genes = item_info)
cat("\n--- [1.5 DGEList 容器檢視 (前 3 樣本文庫大小)] ---\n")
```

    ## 
    ## --- [1.5 DGEList 容器檢視 (前 3 樣本文庫大小)] ---

``` r
print(head(y$samples, 3))
```

    ##                          group lib.size norm.factors
    ## Campus_Normal_R1 Campus.Normal    17297            1
    ## Campus_Normal_R2 Campus.Normal    16585            1
    ## Campus_Normal_R3 Campus.Normal    16403            1

``` r
# 步驟 2. 冷門品項過濾 (Low Sales Filtering: 剔除極低表現雜訊)
```

``` r
# 商業直覺：菜單上有些冷門養生茶或極端特調，全台一天賣不到 1~2 杯。
# 這些品項只會干擾統計模型，無法評估促銷成效，必須事先剔除。
# 標準：保留在至少 3 個樣本中 CPM > 0.5 的品項。

cpm_values <- cpm(y)
keep <- rowSums(cpm_values > 0.5) >= 3

cat("\n--- [步驟 2: 品項過濾統計] ---\n")
```

    ## 
    ## --- [步驟 2: 品項過濾統計] ---

``` r
cat("原始菜單飲品數:", nrow(y), "\n")
```

    ## 原始菜單飲品數: 150

``` r
cat("保留的主力/潛力飲品數:", sum(keep), "\n")
```

    ## 保留的主力/潛力飲品數: 148

``` r
cat("剔除的冷門/零星飲品數:", sum(!keep), "\n")
```

    ## 剔除的冷門/零星飲品數: 2

``` r
# 檢視被剔除的飲品名稱 (看看是不是冷門肉桂冬瓜露、特製丁香苦茶！)
filtered_out_items <- y$genes$DrinkName[!keep]
cat("被剔除的冷門飲品範例:", head(filtered_out_items, 4), "...\n")
```

    ## 被剔除的冷門飲品範例: 銀耳蓮子羹 老薑地瓜露 ...

``` r
# 更新 DGEList 並重設文庫大小 (keep.lib.sizes = FALSE)
y <- y[keep, , keep.lib.sizes = FALSE]
```

``` r
# 步驟 3. 門市總量偏差標準化 (TMM Normalization)
```

``` r
# 商業直覺：促銷日當天全店總杯數暴增（文庫變大），不能只看絕對杯數。
# 若某款茶促銷日賣 100 杯，平日也賣 100 杯，實際上佔比大幅衰退！
# calcNormFactors 使用 TMM (Trimmed Mean of M-values) 算出校正係數 norm.factors。

y <- calcNormFactors(y)
```

    ## calcNormFactors has been renamed to normLibSizes

``` r
cat("\n--- [步驟 3: TMM 標準化係數] ---\n")
```

    ## 
    ## --- [步驟 3: TMM 標準化係數] ---

``` r
print(y$samples[, c("group", "lib.size", "norm.factors")])
```

    ##                          group lib.size norm.factors
    ## Campus_Normal_R1 Campus.Normal    17297    1.1157472
    ## Campus_Normal_R2 Campus.Normal    16585    1.0690894
    ## Campus_Normal_R3 Campus.Normal    16403    1.0704940
    ## Campus_Promo_R1   Campus.Promo    45945    0.8589870
    ## Campus_Promo_R2   Campus.Promo    45611    0.8715179
    ## Campus_Promo_R3   Campus.Promo    43774    0.8683567
    ## Mall_Normal_R1     Mall.Normal    11332    1.2149253
    ## Mall_Normal_R2     Mall.Normal    10446    1.2210037
    ## Mall_Normal_R3     Mall.Normal    10768    1.1687138
    ## Mall_Promo_R1       Mall.Promo    38612    0.8542224
    ## Mall_Promo_R2       Mall.Promo    35057    0.9159824
    ## Mall_Promo_R3       Mall.Promo    34795    0.8880594

``` r
# 步驟 4. 探索性品管分析 (QC & Exploratory Plotting)
```

``` r
# 4.1 MDS Plot (多維尺度分析圖)：檢驗同組重複樣本是否緊密聚集
pdf("plot_01_MDS.pdf", width = 8, height = 7)
col_palette <- c("Campus.Normal" = "darkgreen", "Campus.Promo" = "red",
                 "Mall.Normal" = "blue", "Mall.Promo" = "purple")
pch_shapes <- c("Campus.Normal" = 15, "Campus.Promo" = 16,
                "Mall.Normal" = 17, "Mall.Promo" = 18)

plotMDS(y, col = col_palette[group], pch = pch_shapes[group], cex = 1.3,
        main = "手搖飲門市銷量 MDS 分群圖 (QC)")
```

    ## Warning in title(...): conversion failure on '手搖飲門市銷量 MDS
    ## 分群圖 (QC)' in 'mbcsToSbcs': for 手 (U+624B)

    ## Warning in title(...): conversion failure on '手搖飲門市銷量 MDS
    ## 分群圖 (QC)' in 'mbcsToSbcs': for 搖 (U+6416)

    ## Warning in title(...): conversion failure on '手搖飲門市銷量 MDS
    ## 分群圖 (QC)' in 'mbcsToSbcs': for 飲 (U+98F2)

    ## Warning in title(...): conversion failure on '手搖飲門市銷量 MDS
    ## 分群圖 (QC)' in 'mbcsToSbcs': for 門 (U+9580)

    ## Warning in title(...): conversion failure on '手搖飲門市銷量 MDS
    ## 分群圖 (QC)' in 'mbcsToSbcs': for 市 (U+5E02)

    ## Warning in title(...): conversion failure on '手搖飲門市銷量 MDS
    ## 分群圖 (QC)' in 'mbcsToSbcs': for 銷 (U+92B7)

    ## Warning in title(...): conversion failure on '手搖飲門市銷量 MDS
    ## 分群圖 (QC)' in 'mbcsToSbcs': for 量 (U+91CF)

    ## Warning in title(...): conversion failure on '手搖飲門市銷量 MDS
    ## 分群圖 (QC)' in 'mbcsToSbcs': for 分 (U+5206)

    ## Warning in title(...): conversion failure on '手搖飲門市銷量 MDS
    ## 分群圖 (QC)' in 'mbcsToSbcs': for 群 (U+7FA4)

    ## Warning in title(...): conversion failure on '手搖飲門市銷量 MDS
    ## 分群圖 (QC)' in 'mbcsToSbcs': for 圖 (U+5716)

``` r
legend("topright", legend = levels(group), col = col_palette, pch = pch_shapes, cex = 0.9)
dev.off()
```

    ## png 
    ##   2

``` r
cat("\n>>> 已輸出 plot_01_MDS.pdf (檢查門市與促銷日是否成功分開)\n")
```

    ## 
    ## >>> 已輸出 plot_01_MDS.pdf (檢查門市與促銷日是否成功分開)

``` r
# 4.2 MD Plot (單一樣本表現分佈對稱性)
pdf("plot_02_MD_sample1.pdf", width = 7, height = 6)
plotMD(y, column = 1, main = "樣本 1 (Campus_Normal_R1) Mean-Difference 品管圖")
```

    ## Warning in title(...): conversion failure on '樣本 1
    ## (Campus_Normal_R1) Mean-Difference 品管圖' in 'mbcsToSbcs': for
    ## 樣 (U+6A23)

    ## Warning in title(...): conversion failure on '樣本 1
    ## (Campus_Normal_R1) Mean-Difference 品管圖' in 'mbcsToSbcs': for
    ## 本 (U+672C)

    ## Warning in title(...): conversion failure on '樣本 1
    ## (Campus_Normal_R1) Mean-Difference 品管圖' in 'mbcsToSbcs': for
    ## 品 (U+54C1)

    ## Warning in title(...): conversion failure on '樣本 1
    ## (Campus_Normal_R1) Mean-Difference 品管圖' in 'mbcsToSbcs': for
    ## 管 (U+7BA1)

    ## Warning in title(...): conversion failure on '樣本 1
    ## (Campus_Normal_R1) Mean-Difference 品管圖' in 'mbcsToSbcs': for
    ## 圖 (U+5716)

``` r
abline(h = 0, col = "red", lty = 2, lwd = 2)
dev.off()
```

    ## png 
    ##   2

``` r
cat(">>> 已輸出 plot_02_MD_sample1.pdf\n")
```

    ## >>> 已輸出 plot_02_MD_sample1.pdf

``` r
# 步驟 5. 構建設計矩陣、估計離散度與模型擬合 (Dispersion & QL Fit)
```

``` r
# 5.1 構建無截距設計矩陣 (Cell Means Model)
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)
cat("\n--- [5.1 實驗設計矩陣 (Design Matrix)] ---\n")
```

    ## 
    ## --- [5.1 實驗設計矩陣 (Design Matrix)] ---

``` r
print(design)
```

    ##    Campus.Normal Campus.Promo Mall.Normal Mall.Promo
    ## 1              1            0           0          0
    ## 2              1            0           0          0
    ## 3              1            0           0          0
    ## 4              0            1           0          0
    ## 5              0            1           0          0
    ## 6              0            1           0          0
    ## 7              0            0           1          0
    ## 8              0            0           1          0
    ## 9              0            0           1          0
    ## 10             0            0           0          1
    ## 11             0            0           0          1
    ## 12             0            0           0          1
    ## attr(,"assign")
    ## [1] 1 1 1 1
    ## attr(,"contrasts")
    ## attr(,"contrasts")$group
    ## [1] "contr.treatment"

``` r
# 5.2 估計負二項離散度 (Negative Binomial Dispersion with Empirical Bayes)
# 商業直覺：用全菜單 100 多款飲料的整體波動規律，去「收縮」單品 3 天紀錄的隨機誤差。
y <- estimateDisp(y, design, robust = TRUE)

pdf("plot_03_BCV.pdf", width = 7, height = 6)
plotBCV(y, main = "飲品銷量生物變異係數圖 (BCV Plot)")
```

    ## Warning in title(...): conversion failure on
    ## '飲品銷量生物變異係數圖 (BCV Plot)' in 'mbcsToSbcs': for 飲
    ## (U+98F2)

    ## Warning in title(...): conversion failure on
    ## '飲品銷量生物變異係數圖 (BCV Plot)' in 'mbcsToSbcs': for 品
    ## (U+54C1)

    ## Warning in title(...): conversion failure on
    ## '飲品銷量生物變異係數圖 (BCV Plot)' in 'mbcsToSbcs': for 銷
    ## (U+92B7)

    ## Warning in title(...): conversion failure on
    ## '飲品銷量生物變異係數圖 (BCV Plot)' in 'mbcsToSbcs': for 量
    ## (U+91CF)

    ## Warning in title(...): conversion failure on
    ## '飲品銷量生物變異係數圖 (BCV Plot)' in 'mbcsToSbcs': for 生
    ## (U+751F)

    ## Warning in title(...): conversion failure on
    ## '飲品銷量生物變異係數圖 (BCV Plot)' in 'mbcsToSbcs': for 物
    ## (U+7269)

    ## Warning in title(...): conversion failure on
    ## '飲品銷量生物變異係數圖 (BCV Plot)' in 'mbcsToSbcs': for 變
    ## (U+8B8A)

    ## Warning in title(...): conversion failure on
    ## '飲品銷量生物變異係數圖 (BCV Plot)' in 'mbcsToSbcs': for 異
    ## (U+7570)

    ## Warning in title(...): conversion failure on
    ## '飲品銷量生物變異係數圖 (BCV Plot)' in 'mbcsToSbcs': for 係
    ## (U+4FC2)

    ## Warning in title(...): conversion failure on
    ## '飲品銷量生物變異係數圖 (BCV Plot)' in 'mbcsToSbcs': for 數
    ## (U+6578)

    ## Warning in title(...): conversion failure on
    ## '飲品銷量生物變異係數圖 (BCV Plot)' in 'mbcsToSbcs': for 圖
    ## (U+5716)

``` r
dev.off()
```

    ## png 
    ##   2

``` r
cat(">>> 已輸出 plot_03_BCV.pdf (觀察銷量規模越大、百分比變異越小的規律)\n")
```

    ## >>> 已輸出 plot_03_BCV.pdf (觀察銷量規模越大、百分比變異越小的規律)

``` r
# 5.3 擬合準似然廣義線性模型 (Quasi-Likelihood GLM)
# 商業直覺：每款飲料有自己的脾氣 (例如冰沙容易受天氣影響，熱茶平穩)，給予專屬 QL 離散度。
fit <- glmQLFit(y, design, robust = TRUE)

pdf("plot_04_QLDisp.pdf", width = 7, height = 6)
plotQLDisp(fit, main = "QL 準似然離散度與經驗貝葉斯收縮圖")
```

    ## Warning in title(...): conversion failure on 'QL
    ## 準似然離散度與經驗貝葉斯收縮圖' in 'mbcsToSbcs': for 準 (U+6E96)

    ## Warning in title(...): conversion failure on 'QL
    ## 準似然離散度與經驗貝葉斯收縮圖' in 'mbcsToSbcs': for 似 (U+4F3C)

    ## Warning in title(...): conversion failure on 'QL
    ## 準似然離散度與經驗貝葉斯收縮圖' in 'mbcsToSbcs': for 然 (U+7136)

    ## Warning in title(...): conversion failure on 'QL
    ## 準似然離散度與經驗貝葉斯收縮圖' in 'mbcsToSbcs': for 離 (U+96E2)

    ## Warning in title(...): conversion failure on 'QL
    ## 準似然離散度與經驗貝葉斯收縮圖' in 'mbcsToSbcs': for 散 (U+6563)

    ## Warning in title(...): conversion failure on 'QL
    ## 準似然離散度與經驗貝葉斯收縮圖' in 'mbcsToSbcs': for 度 (U+5EA6)

    ## Warning in title(...): conversion failure on 'QL
    ## 準似然離散度與經驗貝葉斯收縮圖' in 'mbcsToSbcs': for 與 (U+8207)

    ## Warning in title(...): conversion failure on 'QL
    ## 準似然離散度與經驗貝葉斯收縮圖' in 'mbcsToSbcs': for 經 (U+7D93)

    ## Warning in title(...): conversion failure on 'QL
    ## 準似然離散度與經驗貝葉斯收縮圖' in 'mbcsToSbcs': for 驗 (U+9A57)

    ## Warning in title(...): conversion failure on 'QL
    ## 準似然離散度與經驗貝葉斯收縮圖' in 'mbcsToSbcs': for 貝 (U+8C9D)

    ## Warning in title(...): conversion failure on 'QL
    ## 準似然離散度與經驗貝葉斯收縮圖' in 'mbcsToSbcs': for 葉 (U+8449)

    ## Warning in title(...): conversion failure on 'QL
    ## 準似然離散度與經驗貝葉斯收縮圖' in 'mbcsToSbcs': for 斯 (U+65AF)

    ## Warning in title(...): conversion failure on 'QL
    ## 準似然離散度與經驗貝葉斯收縮圖' in 'mbcsToSbcs': for 收 (U+6536)

    ## Warning in title(...): conversion failure on 'QL
    ## 準似然離散度與經驗貝葉斯收縮圖' in 'mbcsToSbcs': for 縮 (U+7E2E)

    ## Warning in title(...): conversion failure on 'QL
    ## 準似然離散度與經驗貝葉斯收縮圖' in 'mbcsToSbcs': for 圖 (U+5716)

``` r
dev.off()
```

    ## png 
    ##   2

``` r
cat(">>> 已輸出 plot_04_QLDisp.pdf (檢視各品項收縮後的最終波動裕度)\n")
```

    ## >>> 已輸出 plot_04_QLDisp.pdf (檢視各品項收縮後的最終波動裕度)

``` r
cat("\n--- [5.3 先驗自由度 (Prior df)] ---\n")
```

    ## 
    ## --- [5.3 先驗自由度 (Prior df)] ---

``` r
print(summary(fit$df.prior))
```

    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    ##   17.67   17.67   17.67   17.67   17.67   17.67

``` r
# 步驟 6. 定義商業對比矩陣與統計檢定 (Contrasts & Testing)
```

``` r
# 6.1 設定四種經典商業考題 (Contrasts)：
#   Q1: 學校門市促銷有沒有爆發？ (Campus.Promo vs Campus.Normal)
#   Q2: 百貨門市促銷有沒有爆發？ (Mall.Promo vs Mall.Normal)
#   Q3: 平常哪間店熱愛買純茶/奶茶？ (Campus.Normal vs Mall.Normal)
#   Q4: 促銷在百貨的威力是否顯著超越學校？ (Interaction / 雙重差異交互作用)
con <- makeContrasts(
  Campus_Promo_Effect = Campus.Promo - Campus.Normal,
  Mall_Promo_Effect   = Mall.Promo - Mall.Normal,
  Store_Baseline_Diff = Campus.Normal - Mall.Normal,
  Interaction_Promo   = (Mall.Promo - Mall.Normal) - (Campus.Promo - Campus.Normal),
  levels = design
)
cat("\n--- [6.1 商業假說對比矩陣] ---\n")
```

    ## 
    ## --- [6.1 商業假說對比矩陣] ---

``` r
print(con)
```

    ##                Contrasts
    ## Levels          Campus_Promo_Effect Mall_Promo_Effect
    ##   Campus.Normal                  -1                 0
    ##   Campus.Promo                    1                 0
    ##   Mall.Normal                     0                -1
    ##   Mall.Promo                      0                 1
    ##                Contrasts
    ## Levels          Store_Baseline_Diff Interaction_Promo
    ##   Campus.Normal                   1                 1
    ##   Campus.Promo                    0                -1
    ##   Mall.Normal                    -1                -1
    ##   Mall.Promo                      0                 1

``` r
# 6.2 傳統 QL F-test 檢定 (檢驗促銷效應是否 != 0)
res_campus <- glmQLFTest(fit, contrast = con[, "Campus_Promo_Effect"])
cat("\n--- [6.2 學校促銷 QL F-test 前 5 大顯著品項] ---\n")
```

    ## 
    ## --- [6.2 學校促銷 QL F-test 前 5 大顯著品項] ---

``` r
print(topTags(res_campus, n = 5))
```

    ## Coefficient:  -1*Campus.Normal 1*Campus.Promo 
    ##              ItemID    DrinkName Category BaseTea BasePrice
    ## DRINK_029 DRINK_029     熟成歐蕾 濃醇奶茶  普洱茶        40
    ## DRINK_035 DRINK_035 仙草凍鮮奶茶 濃醇奶茶  普洱茶        40
    ## DRINK_097 DRINK_097 黑鑽可可冰沙 季節冰沙    紅茶        60
    ## DRINK_028 DRINK_028 黑糖雙波鮮奶 濃醇奶茶  烏龍茶        75
    ## DRINK_027 DRINK_027 經典珍珠奶茶 濃醇奶茶    青茶        35
    ##              logFC   logCPM        F       PValue          FDR
    ## DRINK_029 2.007788 13.74843 72.98855 5.602409e-09 6.154436e-07
    ## DRINK_035 1.858859 13.09164 70.03727 8.316806e-09 6.154436e-07
    ## DRINK_097 1.883667 13.10055 62.89559 2.281178e-08 1.125381e-06
    ## DRINK_028 1.761966 13.89650 59.74372 3.655652e-08 1.316830e-06
    ## DRINK_027 1.759206 13.65842 58.46688 4.448748e-08 1.316830e-06

``` r
# 6.3 總經理的嚴格檢定 glmTreat (門檻: 銷量成長顯著大於 1.5 倍, 即 lfc > log2(1.5))
# 商業直覺：銷量只增加 5% 就算統計顯著也沒商業價值，必須是成長 1.5 倍以上的真爆款！
tr_campus <- glmTreat(fit, contrast = con[, "Campus_Promo_Effect"], lfc = log2(1.5))
cat("\n--- [6.3 glmTreat (成長 > 1.5 倍) 爆款明星飲品 Top 10] ---\n")
```

    ## 
    ## --- [6.3 glmTreat (成長 > 1.5 倍) 爆款明星飲品 Top 10] ---

``` r
top_stars <- topTags(tr_campus, n = 10)
print(top_stars[, c("DrinkName", "Category", "logFC", "PValue", "FDR")])
```

    ## Coefficient:  -1*Campus.Normal 1*Campus.Promo 
    ##                DrinkName Category    logFC       PValue
    ## DRINK_029       熟成歐蕾 濃醇奶茶 2.007788 1.424732e-07
    ## DRINK_035   仙草凍鮮奶茶 濃醇奶茶 1.858859 3.264405e-07
    ## DRINK_097   黑鑽可可冰沙 季節冰沙 1.883667 6.796195e-07
    ## DRINK_106   巨峰紫提雪泥 季節冰沙 1.942657 1.171031e-06
    ## DRINK_028   黑糖雙波鮮奶 濃醇奶茶 1.761966 1.482990e-06
    ## DRINK_027   經典珍珠奶茶 濃醇奶茶 1.759206 1.750380e-06
    ## DRINK_100 奇異果優格冰沙 季節冰沙 1.856845 4.008026e-06
    ## DRINK_102   宇治焙茶碎冰 季節冰沙 1.861080 6.558033e-06
    ## DRINK_099   草莓芝芝雪霜 季節冰沙 1.765854 7.037242e-06
    ## DRINK_040   蜜香燕麥奶茶 濃醇奶茶 1.575077 8.779572e-06
    ##                    FDR
    ## DRINK_029 2.108604e-05
    ## DRINK_035 2.415660e-05
    ## DRINK_097 3.352790e-05
    ## DRINK_106 4.317603e-05
    ## DRINK_028 4.317603e-05
    ## DRINK_027 4.317603e-05
    ## DRINK_100 8.474112e-05
    ## DRINK_102 1.157235e-04
    ## DRINK_099 1.157235e-04
    ## DRINK_040 1.299377e-04

``` r
# 6.4 標記品項狀態 (1: 促銷暴賣, -1: 促銷被排擠滯銷, 0: 無顯著變化)
is_de <- decideTests(tr_campus)
cat("\n--- [6.4 飲品分類統計表] ---\n")
```

    ## 
    ## --- [6.4 飲品分類統計表] ---

``` r
print(summary(is_de))
```

    ##        -1*Campus.Normal 1*Campus.Promo
    ## Down                                 8
    ## NotSig                             113
    ## Up                                  27

``` r
# 6.5 繪製標記差異品項的 MD 圖
pdf("plot_05_MD_DEGs.pdf", width = 8, height = 6)
plotMD(tr_campus, status = is_de, values = c(1, -1), col = c("red", "blue"),
       main = "學校促銷成效 MD 圖 (紅色: 暴賣 > 1.5倍, 藍色: 銷量被排擠)")
```

    ## Warning in title(...): conversion failure on '學校促銷成效 MD 圖
    ## (紅色: 暴賣 > 1.5倍, 藍色: 銷量被排擠)' in 'mbcsToSbcs': for 學
    ## (U+5B78)

    ## Warning in title(...): conversion failure on '學校促銷成效 MD 圖
    ## (紅色: 暴賣 > 1.5倍, 藍色: 銷量被排擠)' in 'mbcsToSbcs': for 校
    ## (U+6821)

    ## Warning in title(...): conversion failure on '學校促銷成效 MD 圖
    ## (紅色: 暴賣 > 1.5倍, 藍色: 銷量被排擠)' in 'mbcsToSbcs': for 促
    ## (U+4FC3)

    ## Warning in title(...): conversion failure on '學校促銷成效 MD 圖
    ## (紅色: 暴賣 > 1.5倍, 藍色: 銷量被排擠)' in 'mbcsToSbcs': for 銷
    ## (U+92B7)

    ## Warning in title(...): conversion failure on '學校促銷成效 MD 圖
    ## (紅色: 暴賣 > 1.5倍, 藍色: 銷量被排擠)' in 'mbcsToSbcs': for 成
    ## (U+6210)

    ## Warning in title(...): conversion failure on '學校促銷成效 MD 圖
    ## (紅色: 暴賣 > 1.5倍, 藍色: 銷量被排擠)' in 'mbcsToSbcs': for 效
    ## (U+6548)

    ## Warning in title(...): conversion failure on '學校促銷成效 MD 圖
    ## (紅色: 暴賣 > 1.5倍, 藍色: 銷量被排擠)' in 'mbcsToSbcs': for 圖
    ## (U+5716)

    ## Warning in title(...): conversion failure on '學校促銷成效 MD 圖
    ## (紅色: 暴賣 > 1.5倍, 藍色: 銷量被排擠)' in 'mbcsToSbcs': for 紅
    ## (U+7D05)

    ## Warning in title(...): conversion failure on '學校促銷成效 MD 圖
    ## (紅色: 暴賣 > 1.5倍, 藍色: 銷量被排擠)' in 'mbcsToSbcs': for 色
    ## (U+8272)

    ## Warning in title(...): conversion failure on '學校促銷成效 MD 圖
    ## (紅色: 暴賣 > 1.5倍, 藍色: 銷量被排擠)' in 'mbcsToSbcs': for 暴
    ## (U+66B4)

    ## Warning in title(...): conversion failure on '學校促銷成效 MD 圖
    ## (紅色: 暴賣 > 1.5倍, 藍色: 銷量被排擠)' in 'mbcsToSbcs': for 賣
    ## (U+8CE3)

    ## Warning in title(...): conversion failure on '學校促銷成效 MD 圖
    ## (紅色: 暴賣 > 1.5倍, 藍色: 銷量被排擠)' in 'mbcsToSbcs': for 倍
    ## (U+500D)

    ## Warning in title(...): conversion failure on '學校促銷成效 MD 圖
    ## (紅色: 暴賣 > 1.5倍, 藍色: 銷量被排擠)' in 'mbcsToSbcs': for 藍
    ## (U+85CD)

    ## Warning in title(...): conversion failure on '學校促銷成效 MD 圖
    ## (紅色: 暴賣 > 1.5倍, 藍色: 銷量被排擠)' in 'mbcsToSbcs': for 色
    ## (U+8272)

    ## Warning in title(...): conversion failure on '學校促銷成效 MD 圖
    ## (紅色: 暴賣 > 1.5倍, 藍色: 銷量被排擠)' in 'mbcsToSbcs': for 銷
    ## (U+92B7)

    ## Warning in title(...): conversion failure on '學校促銷成效 MD 圖
    ## (紅色: 暴賣 > 1.5倍, 藍色: 銷量被排擠)' in 'mbcsToSbcs': for 量
    ## (U+91CF)

    ## Warning in title(...): conversion failure on '學校促銷成效 MD 圖
    ## (紅色: 暴賣 > 1.5倍, 藍色: 銷量被排擠)' in 'mbcsToSbcs': for 被
    ## (U+88AB)

    ## Warning in title(...): conversion failure on '學校促銷成效 MD 圖
    ## (紅色: 暴賣 > 1.5倍, 藍色: 銷量被排擠)' in 'mbcsToSbcs': for 排
    ## (U+6392)

    ## Warning in title(...): conversion failure on '學校促銷成效 MD 圖
    ## (紅色: 暴賣 > 1.5倍, 藍色: 銷量被排擠)' in 'mbcsToSbcs': for 擠
    ## (U+64E0)

``` r
legend("topright", legend = c("顯著暴增 (>1.5x)", "顯著衰退", "無顯著影響"),
       col = c("red", "blue", "black"), pch = 16, cex = 0.8)
```

    ## Warning in text.default(x, y, ...): conversion failure on
    ## '顯著暴增 (>1.5x)' in 'mbcsToSbcs': for 顯 (U+986F)

    ## Warning in text.default(x, y, ...): conversion failure on
    ## '顯著暴增 (>1.5x)' in 'mbcsToSbcs': for 著 (U+8457)

    ## Warning in text.default(x, y, ...): conversion failure on
    ## '顯著暴增 (>1.5x)' in 'mbcsToSbcs': for 暴 (U+66B4)

    ## Warning in text.default(x, y, ...): conversion failure on
    ## '顯著暴增 (>1.5x)' in 'mbcsToSbcs': for 增 (U+589E)

    ## Warning in text.default(x, y, ...): conversion failure on
    ## '顯著衰退' in 'mbcsToSbcs': for 顯 (U+986F)

    ## Warning in text.default(x, y, ...): conversion failure on
    ## '顯著衰退' in 'mbcsToSbcs': for 著 (U+8457)

    ## Warning in text.default(x, y, ...): conversion failure on
    ## '顯著衰退' in 'mbcsToSbcs': for 衰 (U+8870)

    ## Warning in text.default(x, y, ...): conversion failure on
    ## '顯著衰退' in 'mbcsToSbcs': for 退 (U+9000)

    ## Warning in text.default(x, y, ...): conversion failure on
    ## '無顯著影響' in 'mbcsToSbcs': for 無 (U+7121)

    ## Warning in text.default(x, y, ...): conversion failure on
    ## '無顯著影響' in 'mbcsToSbcs': for 顯 (U+986F)

    ## Warning in text.default(x, y, ...): conversion failure on
    ## '無顯著影響' in 'mbcsToSbcs': for 著 (U+8457)

    ## Warning in text.default(x, y, ...): conversion failure on
    ## '無顯著影響' in 'mbcsToSbcs': for 影 (U+5F71)

    ## Warning in text.default(x, y, ...): conversion failure on
    ## '無顯著影響' in 'mbcsToSbcs': for 響 (U+97FF)

``` r
dev.off()
```

    ## png 
    ##   2

``` r
cat(">>> 已輸出 plot_05_MD_DEGs.pdf\n")
```

    ## >>> 已輸出 plot_05_MD_DEGs.pdf

``` r
# 步驟 7. Top 爆款飲品階層分群熱圖 (Heatmap of Sales Winners)
```

``` r
# 計算對數 CPM (log2-CPM, 加入 prior.count=2 避免 log(0))
logCPM <- cpm(y, prior.count = 2, log = TRUE)
rownames(logCPM) <- y$genes$DrinkName

# 挑選在促銷中最顯著變化的前 25 款飲料
top25_indices <- order(tr_campus$table$PValue)[1:25]
top25_mat <- logCPM[top25_indices, ]

# 針對每款飲品進行 Z-score 標準化 (使銷量模式高低對比分明)
top25_scaled <- t(scale(t(top25_mat)))

pdf("plot_06_Heatmap_Winners.pdf", width = 10, height = 8)
col_heatmap <- colorpanel(100, "navy", "white", "firebrick3")
heatmap.2(top25_scaled, col = col_heatmap, scale = "none", trace = "none",
          dendrogram = "both", Rowv = TRUE, Colv = TRUE,
          cexRow = 0.85, cexCol = 1.0, margins = c(10, 12),
          main = "Top 25 促銷反應明星飲品分群熱圖")
```

    ## Warning in axis(4, iy, labels = labRow, las = 2, line = -0.5 +
    ## offsetRow, : conversion failure on '黑鑽可可冰沙' in
    ## 'mbcsToSbcs': for 黑 (U+9ED1)

    ## Warning in axis(4, iy, labels = labRow, las = 2, line = -0.5 +
    ## offsetRow, : conversion failure on '黑鑽可可冰沙' in
    ## 'mbcsToSbcs': for 鑽 (U+947D)

    ## Warning in axis(4, iy, labels = labRow, las = 2, line = -0.5 +
    ## offsetRow, : conversion failure on '黑鑽可可冰沙' in
    ## 'mbcsToSbcs': for 可 (U+53EF)
    ## Warning in axis(4, iy, labels = labRow, las = 2, line = -0.5 +
    ## offsetRow, : conversion failure on '黑鑽可可冰沙' in
    ## 'mbcsToSbcs': for 可 (U+53EF)

    ## Warning in axis(4, iy, labels = labRow, las = 2, line = -0.5 +
    ## offsetRow, : conversion failure on '黑鑽可可冰沙' in
    ## 'mbcsToSbcs': for 冰 (U+51B0)

    ## Warning in axis(4, iy, labels = labRow, las = 2, line = -0.5 +
    ## offsetRow, : conversion failure on '黑鑽可可冰沙' in
    ## 'mbcsToSbcs': for 沙 (U+6C99)

    ## Warning in axis(4, iy, labels = labRow, las = 2, line = -0.5 +
    ## offsetRow, : conversion failure on '奇異果優格冰沙' in
    ## 'mbcsToSbcs': for 奇 (U+5947)

    ## Warning in axis(4, iy, labels = labRow, las = 2, line = -0.5 +
    ## offsetRow, : conversion failure on '奇異果優格冰沙' in
    ## 'mbcsToSbcs': for 異 (U+7570)

    ## Warning in axis(4, iy, labels = labRow, las = 2, line = -0.5 +
    ## offsetRow, : conversion failure on '奇異果優格冰沙' in
    ## 'mbcsToSbcs': for 果 (U+679C)

    ## Warning in axis(4, iy, labels = labRow, las = 2, line = -0.5 +
    ## offsetRow, : conversion failure on '奇異果優格冰沙' in
    ## 'mbcsToSbcs': for 優 (U+512A)

    ## Warning in axis(4, iy, labels = labRow, las = 2, line = -0.5 +
    ## offsetRow, : conversion failure on '奇異果優格冰沙' in
    ## 'mbcsToSbcs': for 格 (U+683C)

    ## Warning in axis(4, iy, labels = labRow, las = 2, line = -0.5 +
    ## offsetRow, : conversion failure on '奇異果優格冰沙' in
    ## 'mbcsToSbcs': for 冰 (U+51B0)

    ## Warning in axis(4, iy, labels = labRow, las = 2, line = -0.5 +
    ## offsetRow, : conversion failure on '奇異果優格冰沙' in
    ## 'mbcsToSbcs': for 沙 (U+6C99)

    ## Warning in axis(4, iy, labels = labRow, las = 2, line = -0.5 +
    ## offsetRow, : conversion failure on '熟成歐蕾' in 'mbcsToSbcs':
    ## for 熟 (U+719F)

    ## Warning in axis(4, iy, labels = labRow, las = 2, line = -0.5 +
    ## offsetRow, : conversion failure on '熟成歐蕾' in 'mbcsToSbcs':
    ## for 成 (U+6210)

    ## Warning in axis(4, iy, labels = labRow, las = 2, line = -0.5 +
    ## offsetRow, : conversion failure on '熟成歐蕾' in 'mbcsToSbcs':
    ## for 歐 (U+6B50)

    ## Warning in axis(4, iy, labels = labRow, las = 2, line = -0.5 +
    ## offsetRow, : conversion failure on '熟成歐蕾' in 'mbcsToSbcs':
    ## for 蕾 (U+857E)

    ## Warning in axis(4, iy, labels = labRow, las = 2, line = -0.5 +
    ## offsetRow, : conversion failure on '伯爵鮮奶茶' in 'mbcsToSbcs':
    ## for 伯 (U+4F2F)

    ## Warning in axis(4, iy, labels = labRow, las = 2, line = -0.5 +
    ## offsetRow, : conversion failure on '伯爵鮮奶茶' in 'mbcsToSbcs':
    ## for 爵 (U+7235)

    ## Warning in axis(4, iy, labels = labRow, las = 2, line = -0.5 +
    ## offsetRow, : conversion failure on '伯爵鮮奶茶' in 'mbcsToSbcs':
    ## for 鮮 (U+9BAE)

    ## Warning in axis(4, iy, labels = labRow, las = 2, line = -0.5 +
    ## offsetRow, : conversion failure on '伯爵鮮奶茶' in 'mbcsToSbcs':
    ## for 奶 (U+5976)

    ## Warning in axis(4, iy, labels = labRow, las = 2, line = -0.5 +
    ## offsetRow, : conversion failure on '伯爵鮮奶茶' in 'mbcsToSbcs':
    ## for 茶 (U+8336)

    ## Warning in axis(4, iy, labels = labRow, las = 2, line = -0.5 +
    ## offsetRow, : conversion failure on '金萱雙Q奶茶' in
    ## 'mbcsToSbcs': for 金 (U+91D1)

    ## Warning in axis(4, iy, labels = labRow, las = 2, line = -0.5 +
    ## offsetRow, : conversion failure on '金萱雙Q奶茶' in
    ## 'mbcsToSbcs': for 萱 (U+8431)

    ## Warning in axis(4, iy, labels = labRow, las = 2, line = -0.5 +
    ## offsetRow, : conversion failure on '金萱雙Q奶茶' in
    ## 'mbcsToSbcs': for 雙 (U+96D9)

    ## Warning in axis(4, iy, labels = labRow, las = 2, line = -0.5 +
    ## offsetRow, : conversion failure on '金萱雙Q奶茶' in
    ## 'mbcsToSbcs': for 奶 (U+5976)

    ## Warning in axis(4, iy, labels = labRow, las = 2, line = -0.5 +
    ## offsetRow, : conversion failure on '金萱雙Q奶茶' in
    ## 'mbcsToSbcs': for 茶 (U+8336)

    ## Warning in axis(4, iy, labels = labRow, las = 2, line = -0.5 +
    ## offsetRow, : conversion failure on '芋圓烏龍奶茶' in
    ## 'mbcsToSbcs': for 芋 (U+828B)

    ## Warning in axis(4, iy, labels = labRow, las = 2, line = -0.5 +
    ## offsetRow, : conversion failure on '芋圓烏龍奶茶' in
    ## 'mbcsToSbcs': for 圓 (U+5713)

    ## Warning in axis(4, iy, labels = labRow, las = 2, line = -0.5 +
    ## offsetRow, : conversion failure on '芋圓烏龍奶茶' in
    ## 'mbcsToSbcs': for 烏 (U+70CF)

    ## Warning in axis(4, iy, labels = labRow, las = 2, line = -0.5 +
    ## offsetRow, : conversion failure on '芋圓烏龍奶茶' in
    ## 'mbcsToSbcs': for 龍 (U+9F8D)

    ## Warning in axis(4, iy, labels = labRow, las = 2, line = -0.5 +
    ## offsetRow, : conversion failure on '芋圓烏龍奶茶' in
    ## 'mbcsToSbcs': for 奶 (U+5976)

    ## Warning in axis(4, iy, labels = labRow, las = 2, line = -0.5 +
    ## offsetRow, : conversion failure on '芋圓烏龍奶茶' in
    ## 'mbcsToSbcs': for 茶 (U+8336)

    ## Warning in axis(4, iy, labels = labRow, las = 2, line = -0.5 +
    ## offsetRow, : conversion failure on '經典珍珠奶茶' in
    ## 'mbcsToSbcs': for 經 (U+7D93)

    ## Warning in axis(4, iy, labels = labRow, las = 2, line = -0.5 +
    ## offsetRow, : conversion failure on '經典珍珠奶茶' in
    ## 'mbcsToSbcs': for 典 (U+5178)

    ## Warning in axis(4, iy, labels = labRow, las = 2, line = -0.5 +
    ## offsetRow, : conversion failure on '經典珍珠奶茶' in
    ## 'mbcsToSbcs': for 珍 (U+73CD)

    ## Warning in axis(4, iy, labels = labRow, las = 2, line = -0.5 +
    ## offsetRow, : conversion failure on '經典珍珠奶茶' in
    ## 'mbcsToSbcs': for 珠 (U+73E0)

    ## Warning in axis(4, iy, labels = labRow, las = 2, line = -0.5 +
    ## offsetRow, : conversion failure on '經典珍珠奶茶' in
    ## 'mbcsToSbcs': for 奶 (U+5976)

    ## Warning in axis(4, iy, labels = labRow, las = 2, line = -0.5 +
    ## offsetRow, : conversion failure on '經典珍珠奶茶' in
    ## 'mbcsToSbcs': for 茶 (U+8336)

    ## Warning in axis(4, iy, labels = labRow, las = 2, line = -0.5 +
    ## offsetRow, : conversion failure on '仙草凍鮮奶茶' in
    ## 'mbcsToSbcs': for 仙 (U+4ED9)

    ## Warning in axis(4, iy, labels = labRow, las = 2, line = -0.5 +
    ## offsetRow, : conversion failure on '仙草凍鮮奶茶' in
    ## 'mbcsToSbcs': for 草 (U+8349)

    ## Warning in axis(4, iy, labels = labRow, las = 2, line = -0.5 +
    ## offsetRow, : conversion failure on '仙草凍鮮奶茶' in
    ## 'mbcsToSbcs': for 凍 (U+51CD)

    ## Warning in axis(4, iy, labels = labRow, las = 2, line = -0.5 +
    ## offsetRow, : conversion failure on '仙草凍鮮奶茶' in
    ## 'mbcsToSbcs': for 鮮 (U+9BAE)

    ## Warning in axis(4, iy, labels = labRow, las = 2, line = -0.5 +
    ## offsetRow, : conversion failure on '仙草凍鮮奶茶' in
    ## 'mbcsToSbcs': for 奶 (U+5976)

    ## Warning in axis(4, iy, labels = labRow, las = 2, line = -0.5 +
    ## offsetRow, : conversion failure on '仙草凍鮮奶茶' in
    ## 'mbcsToSbcs': for 茶 (U+8336)

    ## Warning in axis(4, iy, labels = labRow, las = 2, line = -0.5 +
    ## offsetRow, : conversion failure on '黑糖雙波鮮奶' in
    ## 'mbcsToSbcs': for 黑 (U+9ED1)

    ## Warning in axis(4, iy, labels = labRow, las = 2, line = -0.5 +
    ## offsetRow, : conversion failure on '黑糖雙波鮮奶' in
    ## 'mbcsToSbcs': for 糖 (U+7CD6)

    ## Warning in axis(4, iy, labels = labRow, las = 2, line = -0.5 +
    ## offsetRow, : conversion failure on '黑糖雙波鮮奶' in
    ## 'mbcsToSbcs': for 雙 (U+96D9)

    ## Warning in axis(4, iy, labels = labRow, las = 2, line = -0.5 +
    ## offsetRow, : conversion failure on '黑糖雙波鮮奶' in
    ## 'mbcsToSbcs': for 波 (U+6CE2)

    ## Warning in axis(4, iy, labels = labRow, las = 2, line = -0.5 +
    ## offsetRow, : conversion failure on '黑糖雙波鮮奶' in
    ## 'mbcsToSbcs': for 鮮 (U+9BAE)

    ## Warning in axis(4, iy, labels = labRow, las = 2, line = -0.5 +
    ## offsetRow, : conversion failure on '黑糖雙波鮮奶' in
    ## 'mbcsToSbcs': for 奶 (U+5976)

    ## Warning in axis(4, iy, labels = labRow, las = 2, line = -0.5 +
    ## offsetRow, : conversion failure on '胚芽烤奶茶' in 'mbcsToSbcs':
    ## for 胚 (U+80DA)

    ## Warning in axis(4, iy, labels = labRow, las = 2, line = -0.5 +
    ## offsetRow, : conversion failure on '胚芽烤奶茶' in 'mbcsToSbcs':
    ## for 芽 (U+82BD)

    ## Warning in axis(4, iy, labels = labRow, las = 2, line = -0.5 +
    ## offsetRow, : conversion failure on '胚芽烤奶茶' in 'mbcsToSbcs':
    ## for 烤 (U+70E4)

    ## Warning in axis(4, iy, labels = labRow, las = 2, line = -0.5 +
    ## offsetRow, : conversion failure on '胚芽烤奶茶' in 'mbcsToSbcs':
    ## for 奶 (U+5976)

    ## Warning in axis(4, iy, labels = labRow, las = 2, line = -0.5 +
    ## offsetRow, : conversion failure on '胚芽烤奶茶' in 'mbcsToSbcs':
    ## for 茶 (U+8336)

    ## Warning in axis(4, iy, labels = labRow, las = 2, line = -0.5 +
    ## offsetRow, : conversion failure on '蜜香燕麥奶茶' in
    ## 'mbcsToSbcs': for 蜜 (U+871C)

    ## Warning in axis(4, iy, labels = labRow, las = 2, line = -0.5 +
    ## offsetRow, : conversion failure on '蜜香燕麥奶茶' in
    ## 'mbcsToSbcs': for 香 (U+9999)

    ## Warning in axis(4, iy, labels = labRow, las = 2, line = -0.5 +
    ## offsetRow, : conversion failure on '蜜香燕麥奶茶' in
    ## 'mbcsToSbcs': for 燕 (U+71D5)

    ## Warning in axis(4, iy, labels = labRow, las = 2, line = -0.5 +
    ## offsetRow, : conversion failure on '蜜香燕麥奶茶' in
    ## 'mbcsToSbcs': for 麥 (U+9EA5)

    ## Warning in axis(4, iy, labels = labRow, las = 2, line = -0.5 +
    ## offsetRow, : conversion failure on '蜜香燕麥奶茶' in
    ## 'mbcsToSbcs': for 奶 (U+5976)

    ## Warning in axis(4, iy, labels = labRow, las = 2, line = -0.5 +
    ## offsetRow, : conversion failure on '蜜香燕麥奶茶' in
    ## 'mbcsToSbcs': for 茶 (U+8336)

    ## Warning in axis(4, iy, labels = labRow, las = 2, line = -0.5 +
    ## offsetRow, : conversion failure on '法式可可歐蕾' in
    ## 'mbcsToSbcs': for 法 (U+6CD5)

    ## Warning in axis(4, iy, labels = labRow, las = 2, line = -0.5 +
    ## offsetRow, : conversion failure on '法式可可歐蕾' in
    ## 'mbcsToSbcs': for 式 (U+5F0F)

    ## Warning in axis(4, iy, labels = labRow, las = 2, line = -0.5 +
    ## offsetRow, : conversion failure on '法式可可歐蕾' in
    ## 'mbcsToSbcs': for 可 (U+53EF)
    ## Warning in axis(4, iy, labels = labRow, las = 2, line = -0.5 +
    ## offsetRow, : conversion failure on '法式可可歐蕾' in
    ## 'mbcsToSbcs': for 可 (U+53EF)

    ## Warning in axis(4, iy, labels = labRow, las = 2, line = -0.5 +
    ## offsetRow, : conversion failure on '法式可可歐蕾' in
    ## 'mbcsToSbcs': for 歐 (U+6B50)

    ## Warning in axis(4, iy, labels = labRow, las = 2, line = -0.5 +
    ## offsetRow, : conversion failure on '法式可可歐蕾' in
    ## 'mbcsToSbcs': for 蕾 (U+857E)

    ## Warning in axis(4, iy, labels = labRow, las = 2, line = -0.5 +
    ## offsetRow, : conversion failure on '宇治焙茶碎冰' in
    ## 'mbcsToSbcs': for 宇 (U+5B87)

    ## Warning in axis(4, iy, labels = labRow, las = 2, line = -0.5 +
    ## offsetRow, : conversion failure on '宇治焙茶碎冰' in
    ## 'mbcsToSbcs': for 治 (U+6CBB)

    ## Warning in axis(4, iy, labels = labRow, las = 2, line = -0.5 +
    ## offsetRow, : conversion failure on '宇治焙茶碎冰' in
    ## 'mbcsToSbcs': for 焙 (U+7119)

    ## Warning in axis(4, iy, labels = labRow, las = 2, line = -0.5 +
    ## offsetRow, : conversion failure on '宇治焙茶碎冰' in
    ## 'mbcsToSbcs': for 茶 (U+8336)

    ## Warning in axis(4, iy, labels = labRow, las = 2, line = -0.5 +
    ## offsetRow, : conversion failure on '宇治焙茶碎冰' in
    ## 'mbcsToSbcs': for 碎 (U+788E)

    ## Warning in axis(4, iy, labels = labRow, las = 2, line = -0.5 +
    ## offsetRow, : conversion failure on '宇治焙茶碎冰' in
    ## 'mbcsToSbcs': for 冰 (U+51B0)

    ## Warning in axis(4, iy, labels = labRow, las = 2, line = -0.5 +
    ## offsetRow, : conversion failure on '焦糖布丁冰沙' in
    ## 'mbcsToSbcs': for 焦 (U+7126)

    ## Warning in axis(4, iy, labels = labRow, las = 2, line = -0.5 +
    ## offsetRow, : conversion failure on '焦糖布丁冰沙' in
    ## 'mbcsToSbcs': for 糖 (U+7CD6)

    ## Warning in axis(4, iy, labels = labRow, las = 2, line = -0.5 +
    ## offsetRow, : conversion failure on '焦糖布丁冰沙' in
    ## 'mbcsToSbcs': for 布 (U+5E03)

    ## Warning in axis(4, iy, labels = labRow, las = 2, line = -0.5 +
    ## offsetRow, : conversion failure on '焦糖布丁冰沙' in
    ## 'mbcsToSbcs': for 丁 (U+4E01)

    ## Warning in axis(4, iy, labels = labRow, las = 2, line = -0.5 +
    ## offsetRow, : conversion failure on '焦糖布丁冰沙' in
    ## 'mbcsToSbcs': for 冰 (U+51B0)

    ## Warning in axis(4, iy, labels = labRow, las = 2, line = -0.5 +
    ## offsetRow, : conversion failure on '焦糖布丁冰沙' in
    ## 'mbcsToSbcs': for 沙 (U+6C99)

    ## Warning in axis(4, iy, labels = labRow, las = 2, line = -0.5 +
    ## offsetRow, : conversion failure on '椰香鳳梨冰泥' in
    ## 'mbcsToSbcs': for 椰 (U+6930)

    ## Warning in axis(4, iy, labels = labRow, las = 2, line = -0.5 +
    ## offsetRow, : conversion failure on '椰香鳳梨冰泥' in
    ## 'mbcsToSbcs': for 香 (U+9999)

    ## Warning in axis(4, iy, labels = labRow, las = 2, line = -0.5 +
    ## offsetRow, : conversion failure on '椰香鳳梨冰泥' in
    ## 'mbcsToSbcs': for 鳳 (U+9CF3)

    ## Warning in axis(4, iy, labels = labRow, las = 2, line = -0.5 +
    ## offsetRow, : conversion failure on '椰香鳳梨冰泥' in
    ## 'mbcsToSbcs': for 梨 (U+68A8)

    ## Warning in axis(4, iy, labels = labRow, las = 2, line = -0.5 +
    ## offsetRow, : conversion failure on '椰香鳳梨冰泥' in
    ## 'mbcsToSbcs': for 冰 (U+51B0)

    ## Warning in axis(4, iy, labels = labRow, las = 2, line = -0.5 +
    ## offsetRow, : conversion failure on '椰香鳳梨冰泥' in
    ## 'mbcsToSbcs': for 泥 (U+6CE5)

    ## Warning in axis(4, iy, labels = labRow, las = 2, line = -0.5 +
    ## offsetRow, : conversion failure on '愛文芒果極光冰沙' in
    ## 'mbcsToSbcs': for 愛 (U+611B)

    ## Warning in axis(4, iy, labels = labRow, las = 2, line = -0.5 +
    ## offsetRow, : conversion failure on '愛文芒果極光冰沙' in
    ## 'mbcsToSbcs': for 文 (U+6587)

    ## Warning in axis(4, iy, labels = labRow, las = 2, line = -0.5 +
    ## offsetRow, : conversion failure on '愛文芒果極光冰沙' in
    ## 'mbcsToSbcs': for 芒 (U+8292)

    ## Warning in axis(4, iy, labels = labRow, las = 2, line = -0.5 +
    ## offsetRow, : conversion failure on '愛文芒果極光冰沙' in
    ## 'mbcsToSbcs': for 果 (U+679C)

    ## Warning in axis(4, iy, labels = labRow, las = 2, line = -0.5 +
    ## offsetRow, : conversion failure on '愛文芒果極光冰沙' in
    ## 'mbcsToSbcs': for 極 (U+6975)

    ## Warning in axis(4, iy, labels = labRow, las = 2, line = -0.5 +
    ## offsetRow, : conversion failure on '愛文芒果極光冰沙' in
    ## 'mbcsToSbcs': for 光 (U+5149)

    ## Warning in axis(4, iy, labels = labRow, las = 2, line = -0.5 +
    ## offsetRow, : conversion failure on '愛文芒果極光冰沙' in
    ## 'mbcsToSbcs': for 冰 (U+51B0)

    ## Warning in axis(4, iy, labels = labRow, las = 2, line = -0.5 +
    ## offsetRow, : conversion failure on '愛文芒果極光冰沙' in
    ## 'mbcsToSbcs': for 沙 (U+6C99)

    ## Warning in axis(4, iy, labels = labRow, las = 2, line = -0.5 +
    ## offsetRow, : conversion failure on '巨峰紫提雪泥' in
    ## 'mbcsToSbcs': for 巨 (U+5DE8)

    ## Warning in axis(4, iy, labels = labRow, las = 2, line = -0.5 +
    ## offsetRow, : conversion failure on '巨峰紫提雪泥' in
    ## 'mbcsToSbcs': for 峰 (U+5CF0)

    ## Warning in axis(4, iy, labels = labRow, las = 2, line = -0.5 +
    ## offsetRow, : conversion failure on '巨峰紫提雪泥' in
    ## 'mbcsToSbcs': for 紫 (U+7D2B)

    ## Warning in axis(4, iy, labels = labRow, las = 2, line = -0.5 +
    ## offsetRow, : conversion failure on '巨峰紫提雪泥' in
    ## 'mbcsToSbcs': for 提 (U+63D0)

    ## Warning in axis(4, iy, labels = labRow, las = 2, line = -0.5 +
    ## offsetRow, : conversion failure on '巨峰紫提雪泥' in
    ## 'mbcsToSbcs': for 雪 (U+96EA)

    ## Warning in axis(4, iy, labels = labRow, las = 2, line = -0.5 +
    ## offsetRow, : conversion failure on '巨峰紫提雪泥' in
    ## 'mbcsToSbcs': for 泥 (U+6CE5)

    ## Warning in axis(4, iy, labels = labRow, las = 2, line = -0.5 +
    ## offsetRow, : conversion failure on '草莓芝芝雪霜' in
    ## 'mbcsToSbcs': for 草 (U+8349)

    ## Warning in axis(4, iy, labels = labRow, las = 2, line = -0.5 +
    ## offsetRow, : conversion failure on '草莓芝芝雪霜' in
    ## 'mbcsToSbcs': for 莓 (U+8393)

    ## Warning in axis(4, iy, labels = labRow, las = 2, line = -0.5 +
    ## offsetRow, : conversion failure on '草莓芝芝雪霜' in
    ## 'mbcsToSbcs': for 芝 (U+829D)
    ## Warning in axis(4, iy, labels = labRow, las = 2, line = -0.5 +
    ## offsetRow, : conversion failure on '草莓芝芝雪霜' in
    ## 'mbcsToSbcs': for 芝 (U+829D)

    ## Warning in axis(4, iy, labels = labRow, las = 2, line = -0.5 +
    ## offsetRow, : conversion failure on '草莓芝芝雪霜' in
    ## 'mbcsToSbcs': for 雪 (U+96EA)

    ## Warning in axis(4, iy, labels = labRow, las = 2, line = -0.5 +
    ## offsetRow, : conversion failure on '草莓芝芝雪霜' in
    ## 'mbcsToSbcs': for 霜 (U+971C)

    ## Warning in axis(4, iy, labels = labRow, las = 2, line = -0.5 +
    ## offsetRow, : conversion failure on '荔枝玫瑰凍冰沙' in
    ## 'mbcsToSbcs': for 荔 (U+8354)

    ## Warning in axis(4, iy, labels = labRow, las = 2, line = -0.5 +
    ## offsetRow, : conversion failure on '荔枝玫瑰凍冰沙' in
    ## 'mbcsToSbcs': for 枝 (U+679D)

    ## Warning in axis(4, iy, labels = labRow, las = 2, line = -0.5 +
    ## offsetRow, : conversion failure on '荔枝玫瑰凍冰沙' in
    ## 'mbcsToSbcs': for 玫 (U+73AB)

    ## Warning in axis(4, iy, labels = labRow, las = 2, line = -0.5 +
    ## offsetRow, : conversion failure on '荔枝玫瑰凍冰沙' in
    ## 'mbcsToSbcs': for 瑰 (U+7470)

    ## Warning in axis(4, iy, labels = labRow, las = 2, line = -0.5 +
    ## offsetRow, : conversion failure on '荔枝玫瑰凍冰沙' in
    ## 'mbcsToSbcs': for 凍 (U+51CD)

    ## Warning in axis(4, iy, labels = labRow, las = 2, line = -0.5 +
    ## offsetRow, : conversion failure on '荔枝玫瑰凍冰沙' in
    ## 'mbcsToSbcs': for 冰 (U+51B0)

    ## Warning in axis(4, iy, labels = labRow, las = 2, line = -0.5 +
    ## offsetRow, : conversion failure on '荔枝玫瑰凍冰沙' in
    ## 'mbcsToSbcs': for 沙 (U+6C99)

    ## Warning in axis(4, iy, labels = labRow, las = 2, line = -0.5 +
    ## offsetRow, : conversion failure on '椰果烏龍奶' in 'mbcsToSbcs':
    ## for 椰 (U+6930)

    ## Warning in axis(4, iy, labels = labRow, las = 2, line = -0.5 +
    ## offsetRow, : conversion failure on '椰果烏龍奶' in 'mbcsToSbcs':
    ## for 果 (U+679C)

    ## Warning in axis(4, iy, labels = labRow, las = 2, line = -0.5 +
    ## offsetRow, : conversion failure on '椰果烏龍奶' in 'mbcsToSbcs':
    ## for 烏 (U+70CF)

    ## Warning in axis(4, iy, labels = labRow, las = 2, line = -0.5 +
    ## offsetRow, : conversion failure on '椰果烏龍奶' in 'mbcsToSbcs':
    ## for 龍 (U+9F8D)

    ## Warning in axis(4, iy, labels = labRow, las = 2, line = -0.5 +
    ## offsetRow, : conversion failure on '椰果烏龍奶' in 'mbcsToSbcs':
    ## for 奶 (U+5976)

    ## Warning in axis(4, iy, labels = labRow, las = 2, line = -0.5 +
    ## offsetRow, : conversion failure on '抹茶紅豆鮮奶' in
    ## 'mbcsToSbcs': for 抹 (U+62B9)

    ## Warning in axis(4, iy, labels = labRow, las = 2, line = -0.5 +
    ## offsetRow, : conversion failure on '抹茶紅豆鮮奶' in
    ## 'mbcsToSbcs': for 茶 (U+8336)

    ## Warning in axis(4, iy, labels = labRow, las = 2, line = -0.5 +
    ## offsetRow, : conversion failure on '抹茶紅豆鮮奶' in
    ## 'mbcsToSbcs': for 紅 (U+7D05)

    ## Warning in axis(4, iy, labels = labRow, las = 2, line = -0.5 +
    ## offsetRow, : conversion failure on '抹茶紅豆鮮奶' in
    ## 'mbcsToSbcs': for 豆 (U+8C46)

    ## Warning in axis(4, iy, labels = labRow, las = 2, line = -0.5 +
    ## offsetRow, : conversion failure on '抹茶紅豆鮮奶' in
    ## 'mbcsToSbcs': for 鮮 (U+9BAE)

    ## Warning in axis(4, iy, labels = labRow, las = 2, line = -0.5 +
    ## offsetRow, : conversion failure on '抹茶紅豆鮮奶' in
    ## 'mbcsToSbcs': for 奶 (U+5976)

    ## Warning in axis(4, iy, labels = labRow, las = 2, line = -0.5 +
    ## offsetRow, : conversion failure on '焦糖黑鑽奶茶' in
    ## 'mbcsToSbcs': for 焦 (U+7126)

    ## Warning in axis(4, iy, labels = labRow, las = 2, line = -0.5 +
    ## offsetRow, : conversion failure on '焦糖黑鑽奶茶' in
    ## 'mbcsToSbcs': for 糖 (U+7CD6)

    ## Warning in axis(4, iy, labels = labRow, las = 2, line = -0.5 +
    ## offsetRow, : conversion failure on '焦糖黑鑽奶茶' in
    ## 'mbcsToSbcs': for 黑 (U+9ED1)

    ## Warning in axis(4, iy, labels = labRow, las = 2, line = -0.5 +
    ## offsetRow, : conversion failure on '焦糖黑鑽奶茶' in
    ## 'mbcsToSbcs': for 鑽 (U+947D)

    ## Warning in axis(4, iy, labels = labRow, las = 2, line = -0.5 +
    ## offsetRow, : conversion failure on '焦糖黑鑽奶茶' in
    ## 'mbcsToSbcs': for 奶 (U+5976)

    ## Warning in axis(4, iy, labels = labRow, las = 2, line = -0.5 +
    ## offsetRow, : conversion failure on '焦糖黑鑽奶茶' in
    ## 'mbcsToSbcs': for 茶 (U+8336)

    ## Warning in axis(4, iy, labels = labRow, las = 2, line = -0.5 +
    ## offsetRow, : conversion failure on '黑糖寒天鮮奶茶' in
    ## 'mbcsToSbcs': for 黑 (U+9ED1)

    ## Warning in axis(4, iy, labels = labRow, las = 2, line = -0.5 +
    ## offsetRow, : conversion failure on '黑糖寒天鮮奶茶' in
    ## 'mbcsToSbcs': for 糖 (U+7CD6)

    ## Warning in axis(4, iy, labels = labRow, las = 2, line = -0.5 +
    ## offsetRow, : conversion failure on '黑糖寒天鮮奶茶' in
    ## 'mbcsToSbcs': for 寒 (U+5BD2)

    ## Warning in axis(4, iy, labels = labRow, las = 2, line = -0.5 +
    ## offsetRow, : conversion failure on '黑糖寒天鮮奶茶' in
    ## 'mbcsToSbcs': for 天 (U+5929)

    ## Warning in axis(4, iy, labels = labRow, las = 2, line = -0.5 +
    ## offsetRow, : conversion failure on '黑糖寒天鮮奶茶' in
    ## 'mbcsToSbcs': for 鮮 (U+9BAE)

    ## Warning in axis(4, iy, labels = labRow, las = 2, line = -0.5 +
    ## offsetRow, : conversion failure on '黑糖寒天鮮奶茶' in
    ## 'mbcsToSbcs': for 奶 (U+5976)

    ## Warning in axis(4, iy, labels = labRow, las = 2, line = -0.5 +
    ## offsetRow, : conversion failure on '黑糖寒天鮮奶茶' in
    ## 'mbcsToSbcs': for 茶 (U+8336)

    ## Warning in axis(4, iy, labels = labRow, las = 2, line = -0.5 +
    ## offsetRow, : conversion failure on '小芋圓絲襪奶' in
    ## 'mbcsToSbcs': for 小 (U+5C0F)

    ## Warning in axis(4, iy, labels = labRow, las = 2, line = -0.5 +
    ## offsetRow, : conversion failure on '小芋圓絲襪奶' in
    ## 'mbcsToSbcs': for 芋 (U+828B)

    ## Warning in axis(4, iy, labels = labRow, las = 2, line = -0.5 +
    ## offsetRow, : conversion failure on '小芋圓絲襪奶' in
    ## 'mbcsToSbcs': for 圓 (U+5713)

    ## Warning in axis(4, iy, labels = labRow, las = 2, line = -0.5 +
    ## offsetRow, : conversion failure on '小芋圓絲襪奶' in
    ## 'mbcsToSbcs': for 絲 (U+7D72)

    ## Warning in axis(4, iy, labels = labRow, las = 2, line = -0.5 +
    ## offsetRow, : conversion failure on '小芋圓絲襪奶' in
    ## 'mbcsToSbcs': for 襪 (U+896A)

    ## Warning in axis(4, iy, labels = labRow, las = 2, line = -0.5 +
    ## offsetRow, : conversion failure on '小芋圓絲襪奶' in
    ## 'mbcsToSbcs': for 奶 (U+5976)

    ## Warning in axis(4, iy, labels = labRow, las = 2, line = -0.5 +
    ## offsetRow, : conversion failure on '香草冰淇淋紅茶' in
    ## 'mbcsToSbcs': for 香 (U+9999)

    ## Warning in axis(4, iy, labels = labRow, las = 2, line = -0.5 +
    ## offsetRow, : conversion failure on '香草冰淇淋紅茶' in
    ## 'mbcsToSbcs': for 草 (U+8349)

    ## Warning in axis(4, iy, labels = labRow, las = 2, line = -0.5 +
    ## offsetRow, : conversion failure on '香草冰淇淋紅茶' in
    ## 'mbcsToSbcs': for 冰 (U+51B0)

    ## Warning in axis(4, iy, labels = labRow, las = 2, line = -0.5 +
    ## offsetRow, : conversion failure on '香草冰淇淋紅茶' in
    ## 'mbcsToSbcs': for 淇 (U+6DC7)

    ## Warning in axis(4, iy, labels = labRow, las = 2, line = -0.5 +
    ## offsetRow, : conversion failure on '香草冰淇淋紅茶' in
    ## 'mbcsToSbcs': for 淋 (U+6DCB)

    ## Warning in axis(4, iy, labels = labRow, las = 2, line = -0.5 +
    ## offsetRow, : conversion failure on '香草冰淇淋紅茶' in
    ## 'mbcsToSbcs': for 紅 (U+7D05)

    ## Warning in axis(4, iy, labels = labRow, las = 2, line = -0.5 +
    ## offsetRow, : conversion failure on '香草冰淇淋紅茶' in
    ## 'mbcsToSbcs': for 茶 (U+8336)

    ## Warning in title(main, cex.main = 1.5 * op[["cex.main"]]):
    ## conversion failure on 'Top 25 促銷反應明星飲品分群熱圖' in
    ## 'mbcsToSbcs': for 促 (U+4FC3)

    ## Warning in title(main, cex.main = 1.5 * op[["cex.main"]]):
    ## conversion failure on 'Top 25 促銷反應明星飲品分群熱圖' in
    ## 'mbcsToSbcs': for 銷 (U+92B7)

    ## Warning in title(main, cex.main = 1.5 * op[["cex.main"]]):
    ## conversion failure on 'Top 25 促銷反應明星飲品分群熱圖' in
    ## 'mbcsToSbcs': for 反 (U+53CD)

    ## Warning in title(main, cex.main = 1.5 * op[["cex.main"]]):
    ## conversion failure on 'Top 25 促銷反應明星飲品分群熱圖' in
    ## 'mbcsToSbcs': for 應 (U+61C9)

    ## Warning in title(main, cex.main = 1.5 * op[["cex.main"]]):
    ## conversion failure on 'Top 25 促銷反應明星飲品分群熱圖' in
    ## 'mbcsToSbcs': for 明 (U+660E)

    ## Warning in title(main, cex.main = 1.5 * op[["cex.main"]]):
    ## conversion failure on 'Top 25 促銷反應明星飲品分群熱圖' in
    ## 'mbcsToSbcs': for 星 (U+661F)

    ## Warning in title(main, cex.main = 1.5 * op[["cex.main"]]):
    ## conversion failure on 'Top 25 促銷反應明星飲品分群熱圖' in
    ## 'mbcsToSbcs': for 飲 (U+98F2)

    ## Warning in title(main, cex.main = 1.5 * op[["cex.main"]]):
    ## conversion failure on 'Top 25 促銷反應明星飲品分群熱圖' in
    ## 'mbcsToSbcs': for 品 (U+54C1)

    ## Warning in title(main, cex.main = 1.5 * op[["cex.main"]]):
    ## conversion failure on 'Top 25 促銷反應明星飲品分群熱圖' in
    ## 'mbcsToSbcs': for 分 (U+5206)

    ## Warning in title(main, cex.main = 1.5 * op[["cex.main"]]):
    ## conversion failure on 'Top 25 促銷反應明星飲品分群熱圖' in
    ## 'mbcsToSbcs': for 群 (U+7FA4)

    ## Warning in title(main, cex.main = 1.5 * op[["cex.main"]]):
    ## conversion failure on 'Top 25 促銷反應明星飲品分群熱圖' in
    ## 'mbcsToSbcs': for 熱 (U+71B1)

    ## Warning in title(main, cex.main = 1.5 * op[["cex.main"]]):
    ## conversion failure on 'Top 25 促銷反應明星飲品分群熱圖' in
    ## 'mbcsToSbcs': for 圖 (U+5716)

``` r
dev.off()
```

    ## png 
    ##   2

``` r
cat(">>> 已輸出 plot_06_Heatmap_Winners.pdf\n")
```

    ## >>> 已輸出 plot_06_Heatmap_Winners.pdf

``` r
# 步驟 8. 品類集合富集檢定 (Category Set Enrichment / CAMERA: 相當於 GSEA)
```

``` r
# 商業直覺：生信中的 GO/KEGG Pathway，在飲料店就是「品類 (Category)」！
# 我們想知道：到底「季節冰沙」是全品類集體暴動，還是「原葉純茶」全面被犧牲？
# camera() 是競爭型檢定 (Competitive Test)，會自動校正品類內部飲品間的相互關聯。

# 建立品類清單索引 (相當於 MSigDB gene sets)
cat_list <- split(1:nrow(y), y$genes$Category)

cat("\n--- [步驟 8: 各飲料品類包含的品項數] ---\n")
```

    ## 
    ## --- [步驟 8: 各飲料品類包含的品項數] ---

``` r
print(sapply(cat_list, length))
```

    ## 季節冰沙 原葉純茶 養生溫飲 濃醇奶茶 鮮果特調 
    ##       30       25       23       35       35

``` r
# 執行 camera 檢定
cam_result <- camera(y, index = cat_list, design = design,
                     contrast = con[, "Campus_Promo_Effect"], inter.gene.cor = 0.01)
cat("\n--- [8.1 飲料品類富集檢定結果 (CAMERA)] ---\n")
```

    ## 
    ## --- [8.1 飲料品類富集檢定結果 (CAMERA)] ---

``` r
print(cam_result)
```

    ##          NGenes Direction       PValue          FDR
    ## 季節冰沙     30        Up 8.833651e-06 4.416826e-05
    ## 鮮果特調     35      Down 7.108635e-03 1.777159e-02
    ## 原葉純茶     25      Down 1.056021e-01 1.760035e-01
    ## 養生溫飲     23      Down 4.763321e-01 5.136799e-01
    ## 濃醇奶茶     35        Up 5.136799e-01 5.136799e-01

``` r
# 8.2 Barcode Plot: 視覺化特定品類在銷量成長光譜上的集中趨勢
pdf("plot_07_Barcode_Smoothies.pdf", width = 8, height = 5)
res_for_barcode <- glmQLFTest(fit, contrast = con[, "Campus_Promo_Effect"])
smoothie_idx <- cat_list[["季節冰沙"]]
barcodeplot(res_for_barcode$table$logFC, index = smoothie_idx,
            labels = c("促銷大跌", "促銷暴增"),
            main = "季節冰沙品類於促銷全體光譜之富集條形圖 (Barcode Plot)")
```

    ## Warning in title(...): conversion failure on
    ## '季節冰沙品類於促銷全體光譜之富集條形圖 (Barcode Plot)' in
    ## 'mbcsToSbcs': for 季 (U+5B63)

    ## Warning in title(...): conversion failure on
    ## '季節冰沙品類於促銷全體光譜之富集條形圖 (Barcode Plot)' in
    ## 'mbcsToSbcs': for 節 (U+7BC0)

    ## Warning in title(...): conversion failure on
    ## '季節冰沙品類於促銷全體光譜之富集條形圖 (Barcode Plot)' in
    ## 'mbcsToSbcs': for 冰 (U+51B0)

    ## Warning in title(...): conversion failure on
    ## '季節冰沙品類於促銷全體光譜之富集條形圖 (Barcode Plot)' in
    ## 'mbcsToSbcs': for 沙 (U+6C99)

    ## Warning in title(...): conversion failure on
    ## '季節冰沙品類於促銷全體光譜之富集條形圖 (Barcode Plot)' in
    ## 'mbcsToSbcs': for 品 (U+54C1)

    ## Warning in title(...): conversion failure on
    ## '季節冰沙品類於促銷全體光譜之富集條形圖 (Barcode Plot)' in
    ## 'mbcsToSbcs': for 類 (U+985E)

    ## Warning in title(...): conversion failure on
    ## '季節冰沙品類於促銷全體光譜之富集條形圖 (Barcode Plot)' in
    ## 'mbcsToSbcs': for 於 (U+65BC)

    ## Warning in title(...): conversion failure on
    ## '季節冰沙品類於促銷全體光譜之富集條形圖 (Barcode Plot)' in
    ## 'mbcsToSbcs': for 促 (U+4FC3)

    ## Warning in title(...): conversion failure on
    ## '季節冰沙品類於促銷全體光譜之富集條形圖 (Barcode Plot)' in
    ## 'mbcsToSbcs': for 銷 (U+92B7)

    ## Warning in title(...): conversion failure on
    ## '季節冰沙品類於促銷全體光譜之富集條形圖 (Barcode Plot)' in
    ## 'mbcsToSbcs': for 全 (U+5168)

    ## Warning in title(...): conversion failure on
    ## '季節冰沙品類於促銷全體光譜之富集條形圖 (Barcode Plot)' in
    ## 'mbcsToSbcs': for 體 (U+9AD4)

    ## Warning in title(...): conversion failure on
    ## '季節冰沙品類於促銷全體光譜之富集條形圖 (Barcode Plot)' in
    ## 'mbcsToSbcs': for 光 (U+5149)

    ## Warning in title(...): conversion failure on
    ## '季節冰沙品類於促銷全體光譜之富集條形圖 (Barcode Plot)' in
    ## 'mbcsToSbcs': for 譜 (U+8B5C)

    ## Warning in title(...): conversion failure on
    ## '季節冰沙品類於促銷全體光譜之富集條形圖 (Barcode Plot)' in
    ## 'mbcsToSbcs': for 之 (U+4E4B)

    ## Warning in title(...): conversion failure on
    ## '季節冰沙品類於促銷全體光譜之富集條形圖 (Barcode Plot)' in
    ## 'mbcsToSbcs': for 富 (U+5BCC)

    ## Warning in title(...): conversion failure on
    ## '季節冰沙品類於促銷全體光譜之富集條形圖 (Barcode Plot)' in
    ## 'mbcsToSbcs': for 集 (U+96C6)

    ## Warning in title(...): conversion failure on
    ## '季節冰沙品類於促銷全體光譜之富集條形圖 (Barcode Plot)' in
    ## 'mbcsToSbcs': for 條 (U+689D)

    ## Warning in title(...): conversion failure on
    ## '季節冰沙品類於促銷全體光譜之富集條形圖 (Barcode Plot)' in
    ## 'mbcsToSbcs': for 形 (U+5F62)

    ## Warning in title(...): conversion failure on
    ## '季節冰沙品類於促銷全體光譜之富集條形圖 (Barcode Plot)' in
    ## 'mbcsToSbcs': for 圖 (U+5716)

    ## Warning in axis(side = 2, at = lab.at, padj = 3.8, cex.axis =
    ## 0.85, labels = labels[1], : conversion failure on '促銷大跌' in
    ## 'mbcsToSbcs': for 促 (U+4FC3)

    ## Warning in axis(side = 2, at = lab.at, padj = 3.8, cex.axis =
    ## 0.85, labels = labels[1], : conversion failure on '促銷大跌' in
    ## 'mbcsToSbcs': for 銷 (U+92B7)

    ## Warning in axis(side = 2, at = lab.at, padj = 3.8, cex.axis =
    ## 0.85, labels = labels[1], : conversion failure on '促銷大跌' in
    ## 'mbcsToSbcs': for 大 (U+5927)

    ## Warning in axis(side = 2, at = lab.at, padj = 3.8, cex.axis =
    ## 0.85, labels = labels[1], : conversion failure on '促銷大跌' in
    ## 'mbcsToSbcs': for 跌 (U+8DCC)

    ## Warning in axis(side = 4, at = lab.at, padj = -3.8, cex.axis =
    ## 0.85, labels = labels[2], : conversion failure on '促銷暴增' in
    ## 'mbcsToSbcs': for 促 (U+4FC3)

    ## Warning in axis(side = 4, at = lab.at, padj = -3.8, cex.axis =
    ## 0.85, labels = labels[2], : conversion failure on '促銷暴增' in
    ## 'mbcsToSbcs': for 銷 (U+92B7)

    ## Warning in axis(side = 4, at = lab.at, padj = -3.8, cex.axis =
    ## 0.85, labels = labels[2], : conversion failure on '促銷暴增' in
    ## 'mbcsToSbcs': for 暴 (U+66B4)

    ## Warning in axis(side = 4, at = lab.at, padj = -3.8, cex.axis =
    ## 0.85, labels = labels[2], : conversion failure on '促銷暴增' in
    ## 'mbcsToSbcs': for 增 (U+589E)

``` r
dev.off()
```

    ## png 
    ##   2

``` r
cat(">>> 已輸出 plot_07_Barcode_Smoothies.pdf\n")
```

    ## >>> 已輸出 plot_07_Barcode_Smoothies.pdf

``` r
# 學生課後實作與思考題 (Student Practice Tasks)
```

``` r
cat("\n================================================================\n")
```

    ## 
    ## ================================================================

``` r
cat("【學生實作作業 Challenge Tasks】\n")
```

    ## 【學生實作作業 Challenge Tasks】

``` r
cat("1. [練習一: 百貨門市促銷分析]\n")
```

    ## 1. [練習一: 百貨門市促銷分析]

``` r
cat("   請使用 tr_mall <- glmTreat(fit, contrast = con[,'Mall_Promo_Effect'], lfc=log2(1.5))\n")
```

    ##    請使用 tr_mall <- glmTreat(fit, contrast = con[,'Mall_Promo_Effect'], lfc=log2(1.5))

``` r
cat("   輸出百貨門市促銷前 5 大明星飲品，比較與學校門市有何不同？\n\n")
```

    ##    輸出百貨門市促銷前 5 大明星飲品，比較與學校門市有何不同？

``` r
cat("2. [練習二: 平日門市客群分析]\n")
```

    ## 2. [練習二: 平日門市客群分析]

``` r
cat("   請檢驗 con[,'Store_Baseline_Diff']，找出平日學校比百貨更受歡迎的飲品。\n\n")
```

    ##    請檢驗 con[,'Store_Baseline_Diff']，找出平日學校比百貨更受歡迎的飲品。

``` r
cat("3. [練習三: 雙重差異交互作用分析]\n")
```

    ## 3. [練習三: 雙重差異交互作用分析]

``` r
cat("   檢驗 con[,'Interaction_Promo']，找出『在百貨促銷成長幅度顯著大於學校促銷』的特殊飲品！\n\n")
```

    ##    檢驗 con[,'Interaction_Promo']，找出『在百貨促銷成長幅度顯著大於學校促銷』的特殊飲品！

``` r
cat("4. [練習四: 門檻敏感度分析]\n")
```

    ## 4. [練習四: 門檻敏感度分析]

``` r
cat("   如果將 glmTreat 的門檻從 1.5 倍 (log2(1.5)) 調升至 2 倍 (log2(2))，\n")
```

    ##    如果將 glmTreat 的門檻從 1.5 倍 (log2(1.5)) 調升至 2 倍 (log2(2))，

``` r
cat("   學校門市符合爆款資格的飲品數量會如何改變？\n")
```

    ##    學校門市符合爆款資格的飲品數量會如何改變？

``` r
cat("================================================================\n\n")
```

    ## ================================================================

``` r
cat(">>> 全部分析腳本執行完畢！請檢視產生的 7 張 PDF 圖檔與 Console 報表。\n")
```

    ## >>> 全部分析腳本執行完畢！請檢視產生的 7 張 PDF 圖檔與 Console 報表。
