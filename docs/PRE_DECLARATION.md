# 预先声明 — 主分析集迁移 (v1 → v2)

**手稿：** JARE-D-26-04883
**声明日期：** 2026-08-09
**声明时点：** 在 Stage 2 任何计算作业提交之前
**执行目录：** `/home/gao/eccDNA/New_run_202060809`

本文件在重算开始前写定并存档。它的作用是把这一轮重跑固定为一次**数据集纠错**，
而不是看到结果之后的选择性调整。

---

## 1. 变更内容

主分析 call set 由 `circlemap_current` 改为 `circlemap_methods`。

| | v1（原投稿/修回稿） | v2（本轮） |
|---|---|---|
| 主 call set | 归档的 Circle-Map BED，未施加声明的过滤 | 满足全部五条过滤的调用 |
| 总调用数 | 16,857,160 | 10,099,219（60.1%） |
| COVID / HC | 15,006,670 / 1,850,490 | 8,830,916 / 1,268,303 |
| 每样本中位保留率 | 1.000（定义） | 0.654 |

## 2. 变更理由

手稿 Methods 段声明：

> High-confidence candidates met all of the following uniform filters based on our
> previous study: (i) ≥3 split reads supporting the junction; (ii) Circle-Map score
> >200; (iii) mean per-base coverage across the circle exceeding the standard
> deviation of its per-base coverage profile; (iv) a ≥0.3 coverage increase at both
> start and end relative to local flanking windows; and (v) <10% uncovered bases
> within the circle. The final call set (BED) included junction coordinates and
> support metrics for downstream analyses.

对归档 BED 逐行施加这五条后，只有 60.1% 的调用通过。即**全部下游分析实际使用的
数据集从未施加过 Methods 声明的过滤**。失效分布：

- 覆盖度三条（iii–v）剔除 2,058,259 条（12.2%）
- 在此基础上 split ≥3 / score >200 再剔除 4,649,709 条（27.7%）

另有两处内部不一致由同一原因导致：

1. 稳健性阶梯（strict、very_strict、artifact_masked、high_support_masked、
   consensus_t10）全部以 methods 为根构建，只有主集 `current` 游离在外。
2. 四个实验验证靶点的支持数本就是在 methods 集上计算的，与主分析不同源。

## 3. 保持不变的内容

以下在整轮重跑中一律不得变动。任何变动都会使 v1/v2 差异不再可归因于 call set。

- **统计方法**：双侧 Wilcoxon 秩和检验、Benjamini–Hochberg 校正、Cliff's δ、
  Fisher 精确检验、HC3 稳健标准误 OLS。
- **显著性阈值**：BH q < 0.05；差异 eccGene 附加 |log2FC| ≥ 1。
- **基因归属定义**：junction（主）、interval、midpoint 三种定义及其判定规则。
- **EA 归一化公式**：按基因长度并按样本内跨基因求和归一。
- **过滤与掩码参数**：五条过滤阈值、split/score 阶梯、±100 bp 断点窗、
  Umap k=100 ≥ 0.90、blacklist/gap/segdup/RepeatMasker 规则。
- **随机种子与重采样次数**：片段长度 seed 20260729 / 1,000 次自举；
  RCA 5,000 次自举；P0 降采样 100 次重复、目标 4,000 与 20,000。
- **验证靶点**：CTNNA2circle、SDK1-circle、TCF7L1-circle 不重选。三者本就在
  methods 集上选出，且湿实验已排期。
- **methods call set 本身**：只按已记录的 SHA-256 校验，不重建。

## 4. 预先接受的后果

以下差异在声明时已经量化，属于预期结果，不构成调整分析规则的理由。

| 项目 | 已测差异 |
|---|---|
| 受检基因数 | 26,642 → 26,087 |
| 差异 eccGene 数 | 3,402 → 3,487；交集 2,306；Jaccard 0.503 |
| Fig 3b top-30 榜单 | 保留 9/30；LOC101929073、RORA-AS2、TRAF4 掉出 |
| 染色质 absolute burden | 方向 30/30，显著性 30/30（稳定） |
| 染色质 relative enrichment | 方向 24/30，显著性 25/30（5 项翻转） |
| A549 H3K27ac / H3K4me3 | 由不显著转为显著（δ −0.19 → −0.37） |
| CD14 H3K9ac / H3K4me3 | 失去显著性 |
| EPM 中位（COVID vs HC） | 3,143 vs 814 → 2,053 vs 552 |
| 协变量校正倍数 | 5.08 (2.74–9.41) → 4.84 (2.62–8.93) |

尚未计算、结果未知的项目：片段长度峰值（v1 为 196/365/571 bp，间距
169/206 bp）、P0 协变量与降采样对照、Figure 5 位点级表、BCL3/PROCR 排序。

## 5. 处置政策

1. 所有结果按算出来的样子报告。不因结果方向调整任何阈值、定义或种子。
2. 若某项 v1 中显著的结果在 v2 中不再显著，从正文删除或降级为不显著，
   不得改用其他检验或其他定义使其恢复显著。
3. 若 v2 出现 v1 中没有的显著结果（已知 A549 H3K27ac/H3K4me3），如实增补，
   并在 Discussion 中说明它出现在过滤后的数据中。
4. v1 全部结果保留归档，不删除。审稿回复中引用过的数字必须可追溯。
5. 正文与三份审稿回复中的每一个数字，由脚本对照 v2 的表逐条核验，不依赖人工通读。

## 6. 向审稿人的披露

本次变更将在修回信中主动说明，不作静默替换。理由：Reviewer 3 Major Comment 2
已在追问 call set 与检出流程的稳健性，主动交代远好于被发现前后两版数字不一致。

`circlemap_current` 将作为一条**宽松敏感性臂**保留在稳健性分析中并明确报告。

---

**签署时点的 Stage 0 校验结果：** 78/78 个 `circlemap_methods` BED 的 SHA-256
与 `circlemap_output_callset_sample_metrics.tsv` 中记录值全部匹配，无缺失、
无多余、无不符。校验记录见 `logs/S0_callset_verification.txt`。

---

## 修订 A1 — 片段长度等量抽样数（2026-08-09，Stage 2 运行中）

**变更：** 等量抽样目标由每样本 4,000 条唯一 eccDNA 改为 2,500 条。

**理由：不是选择，是算术上不可行。** v1 取 4,000 是因为当时所有样本都不低于该值。
施加五条过滤后，最小样本 `ZXS181` 只剩 2,598 条唯一 eccDNA，分析在
`run_equal_count_subsampling` 处直接抛出 `ValueError` 终止：

```
ValueError: Cannot subsample 4000 unique eccDNAs per sample: ZXS181 has 2598
```

2,500 是全部 78 个样本都能提供的最大整数值。该值在**看到任何 v2 结果之前**、
仅根据可行性下界确定，与结果方向无关。

**未变更：** 随机种子（20260729）、自举次数（1,000）、抽样重复次数（100）、
平滑参数、峰值搜索区间、共识判定规则全部保持 v1 设定。

**已知影响：** 抽样数下降会略微加宽等量抽样臂的置信区间。该臂是敏感性检查，
主结论来自全量密度曲线，不受影响。

**记录：** 作业在失败前已经完成主分析路径，输出的候选峰为 **196 / 366 / 571 /
764 bp**（v1 为 196 / 365 / 571 / 760 bp）。Gate 2 据此实质解除——峰位在
methods 集下未发生实质移动。最终数值以补跑完成后的表为准。

---

## 修订 A2 — 验证靶点集合（2026-08-09，Stage 2 运行中）

**依据：** `功能转染实验_eccDNA_20260806.docx`（现行功能验证方案）。

**最终实验验证靶点为三个**，非手稿 v1 中讨论的四个：

| 靶标 | BED 坐标 (hg38) | 长度 | 去向 |
|---|---|---:|---|
| CTNNA2circle | chr2:79,885,063–79,885,444 | 381 bp | 实验验证 |
| SDK1-circle | chr7:3,619,742–3,620,122 | 380 bp | 实验验证 |
| TCF7L1-circle | chr2:85,213,855–85,214,235 | 380 bp | 实验验证 |
| CAB39circle | chr2:230,727,472–230,727,831 | 359 bp | 仅计算审计 |
| chr22 区间 | chr22:42,359,179–42,359,558 | 379 bp | 仅计算审计 |

**处置原则：计算审计五个，实验验证三个。**

CAB39 与 chr22 区间保留在计算审计中，不删除。理由：Reviewer 1 Major Comment 7
原文点名 `CAB39circle`，手稿必须给出它未进入实验验证的依据（两个断点窗均与注释
重复序列重叠，未通过严格断点 mask）。把它们从审计里删掉会毁掉该决定的证据链。

**TCF7L1 此前从未进入稳健性流水线。** v1 的靶点列表硬编码为四个，分布在四个脚本
中。本轮在 `scripts/hpc/v2_patched/` 下建立打补丁的副本，加入 TCF7L1，v1 原脚本
不作修改：

- `audit_validation_target_masks.py` — TARGETS 增加 TCF7L1
- `run_core_robustness_analysis.py` — TARGETS 增加 TCF7L1（逐 call set 支持数）
- `generate_final_report.py` — 三处靶点/基因元组扩充
- `summarize_downstream_robustness.py` — CANDIDATE_GENES 扩充

补丁只扩大受审靶点范围，不改变任何判定规则、阈值或统计方法。

**本轮已获得的 TCF7L1 独立审计结果**（对全部 78 个原始 BED 直接扫描，
`scripts/local/audit_tcf7l1.py`，以 CTNNA2 与 SDK1 作阳性对照并复现其已知数值）：

- COVID 精确检出 9/39，HC 0/39
- 满足五条过滤 8/39（XS75 因 start coverage increase = 0.155 < 0.3 未通过）
- ±10 bp 容差下仍为 9 例，无非精确邻近断点
- 两个断点窗 Umap k=100 = 1.000，RepeatMasker / segdup / blacklist 重叠均为 0，
  **通过严格断点 mask**

**对 S5 的影响：** 手稿 Results 中"Of the four candidate validation targets…"一段、
Figure S7d,e 与 Table S16 全部需按五审计/三验证的结构重写。

---

## 修订 A3 — 授权重新实现三处缺失的上游分析（2026-08-09，Stage 4 进行中）

**变更：** 三处分析的原始代码不在仓库中，经用户明确授权，按 Methods 描述重新实现。

**涉及范围：**

| 产物 | 缺失的上游代码 | 影响的图表 |
|---|---|---|
| `chromosome_distribution.tsv` | 每样本按染色体归一化 | Figure 2b、Table S10 |
| `repeat_class_enrichment.tsv` | 每样本重复类映射比 | Figure 2c、Table S11 |
| `gene_feature_enrichment_scores.tsv` | 每样本基因元件 O/E | Figure 2d、Table S12/S13 |
| 复发精确区间矩阵 | ≥10/39 COVID 且 HC=0 的精确区间 | Figure 3d |
| `A1–A6_*.tsv` | 临床相关分析 | Figure S5、临床相关段落 |

**这是重新实现，不是复现。** 与本轮其余部分不同，这几项没有"用 v1 输入复现已投稿产物"的保真门可用——因为原始代码不存在，无法验证我的实现与原实现是否一致。所有产物必须在图注、表注与修回信中标注为按 Methods 重新实现。

**从 Methods 与既有表的列名可以精确复原的定义：**

- 基因元件富集分数 = `log2(observed_fraction / expected_fraction)`，已在 v1 表上逐行验证通过。
- `observed_fraction` = 该样本落在该元件的 eccDNA 数 / 该样本 eccDNA 总数。
- `expected_fraction` = 背景集落在该元件的数 / 背景集总数。
- 染色体归一化值 = (该染色体 eccDNA 数 / 样本总数 × 100) / 染色体长度(Mb)。
- 重复类映射比 = (落在该类的 eccDNA 占比) / (该类占基因组的比例)。
- eccDNA 归属按**起始坐标**（1 bp）判定，与 junction 定义一致（Methods 第 82 行）。
- 工具与版本：bedtools v2.28.0、RepeatMasker open-4.0.5、hg38。

**已识别的实现分歧点，逐条记录：**

1. **背景集定义无法复原。** v1 的 `background_total_b = 16,717,103` 是跨全部样本与元件的单一常数，与 current 集总数（16,857,160）、有效行数（16,807,187）、canonical 数（16,857,075）均不相等，其过滤规则无从得知。本轮以 **v2 全队列汇总调用**重建背景，理由是：观测值来自 v2，期望值必须同源，否则是拿 v2 的观测比 v1 的期望。此选择使 v2 的 O/E 与 v1 的 O/E **不可直接逐值比较**，只可比较方向与组间差异。
2. **重叠判定的边界约定未知。** bedtools 默认 1 bp 重叠即计入，但 v1 是否使用 `-f` 最小重叠比例、是否 `-u` 去重、是否限制 canonical 染色体，均无记录。本轮采用：canonical 染色体、起始坐标点落入元件即计入、每个 eccDNA 对每个元件至多计一次。
3. **重复类的分类粒度未知。** RepeatMasker 的 `repClass` 有 16 个值与图注列出的类别一致，但 v1 是否合并了子类（如 `DNA?` 与 `DNA`）无记录。本轮不做合并，并输出实际出现的类别清单以供核对。
4. **元件注释来源版本未固定。** Methods 只写 refGene 与 CpG island，未记录下载日期或具体表。本轮使用与 eccGene 分析同一份 `refGene.txt`，CpG island 使用 UCSC 同期表，并记录其校验和。
5. **临床相关的置信区间方法未知。** A2/A3 的 `ci_lo/ci_hi` 未标注是自举还是 Fisher z 变换。本轮采用 Fisher z 变换（Spearman 的标准做法），并在表中新增一列标明方法。
6. **Figure 3c 的 Metascape 词条无法重新获取。** Metascape 无公开 API，且其结果依赖提交时的数据库版本，即使人工重新提交也不保证可比。该面板不在本次重新实现范围内。

**处置：** 每一项产物都写入独立的 `*_REIMPLEMENTED.md` 说明，列出所用定义、参数与上述分歧点。正文与修回信中相应位置须注明"按 Methods 重新实现"。若原始代码后续被找到，应以原实现重算并替换。

---

## 修订 A4 — 重新实现的口径以 Methods 为准（2026-08-09，Stage 4 进行中）

**背景更正。** 修订 A3 写"原始代码不在仓库中"。这在本地仓库范围内属实，但随后在 HPC 项目树上（`/gpfs/data/gao/Covid-eccDNA/annotation/` 与 `Repeat_annotation/`，均不在 `Revise/` 下）找到了基因元件与重复类的原始脚本。染色体分布与临床相关的脚本经全盘搜索仍未找到。

**用户决定：按 Methods 的方法重新实现。** 在原始代码与 Methods 描述不一致处，以 Methods 为准。

**已确认的两处不一致：**

1. **归属坐标。** Methods 第 82 行："we counted unique eccDNA molecules by their
   **start coordinates** and mapped them to seven genomic elements … using bedtools"。
   原始脚本 `eccDNA_element_ratio.sh` 实际执行的是
   `bedtools intersect -a "$BED" -b "${ELEMENT}_merged.bed" -u`，即按 eccDNA
   **完整区间**判定。脚本注释中还保留着"若只按起始位置统计，可把 -a 替换为
   `<(awk '{print $1,$2,$2+1}')`"，说明按起点的做法当时被考虑过但未采用。
   **本轮按 Methods 采用起始坐标。**

2. **O/E 的分母。** Methods 中的公式（手稿 OMML 公式 [1]）为：

   > O/E ratio of elements = (ratio of unique eccDNAs falling in a given element)
   > / (ratio of the genome occupied by that element)

   即分母是**该元件占基因组的比例**。而 v1 的
   `gene_feature_enrichment_scores.tsv` 用的是全队列汇总 eccDNA 背景
   （`background_overlap_e / background_total_b`，跨样本恒定的
   16,717,103）。两者数值明显不同：exon 的期望比例，Methods 口径为 0.02856，
   v1 表为 0.07007。**本轮按 Methods 采用基因组占比**，取自
   `annotation/genome_element_coverage_ratio.txt`，并由合并后的元件 BED 独立重算校验。

**后果，须在正文与修回信中明写：** v2 的基因元件 O/E 与 v1 的数值**不可直接逐值比较**，
差异同时来自 call set 迁移与口径更正两个来源。方向性结论（哪些元件富集、
组间差异的符号）仍可比较。这一处使 Methods 的描述与实际计算首次一致。

**Figure 2c（重复类）不在重新实现范围内。** 原始脚本 `Repeat_ratio_{covid,normal}.sh`
从 BAM 经 `samtools view | bedtools bamtobed | bedtools intersect` 计算读段级
映射比，完全不经过 eccDNA 调用集。BAM 在 v1 与 v2 之间相同，因此该分析的
全部数值不受本次迁移影响，正文相应段落无需改动。

**Figure 2b（染色体分布）口径：** 按 Methods "normalizing by chromosomal length
(percentage of eccDNAs per Mb)"，即
`(该染色体 eccDNA 数 / 样本 eccDNA 总数 × 100) / 染色体长度(Mb)`，
与 v1 表列 `normalized_fraction_per_mb` 的量纲一致。

**Figure S5（临床相关）** 的原始代码仍未找到，按 A3 第 5 条以 Fisher z 变换
计算 Spearman 置信区间，并在表中标明方法。

**A4 补充（用户确认）：** 基因元件采用**起始坐标**计数，即修订 A4 第 1 条的 Methods 口径。
已知后果：Figure 2d 的 O/E 数值较已发表版本大幅下降（NS120 exon 由 2.94 降至 1.37），
原因是计数约定而非 call set 迁移——完整区间会把一条 eccDNA 计入它跨越的每一个元件，
起始坐标只计一次。正文与修回信须写明这是口径更正。

**独立于本次迁移的正文错误：** Results 中"observed/expected ratios above 1 at all seven
elements, lowest at introns and highest at 5′UTR and exons"一句，在 v2（CpG islands
在 COVID 组为 0.927）与 **v1 自身的表**（CpG islands 0.850）下均不成立，且两版中最低
的都是 CpG islands 而非 intron。该句须独立更正。

---

## 修订 A5 — 算术与三靶点时间线澄清（2026-08-12，投稿一致性复核）

本节是对冻结记录的追加更正，不回写或删除上文的历史内容。

1. **保留比例分母。** `10,099,219 / 16,857,160 = 59.91056%`，按一位小数
   应报告为 **59.9%**。上文的 60.1% 来自以 16,807,187 条有效解析记录为分母
   （60.08869%）。投稿正文与回复若同时引用 16,857,160 的归档总数，应统一写
   59.9%；不得把两个分母混在同一句中。
2. **每样本中位保留率。** 最终锁定表的中位数为 **0.652**；上文 0.654 是重跑前
   的初步汇总值，后续文件以 0.652 为准。
3. **三靶点时间线。** `功能转染实验_eccDNA_20260806.docx` 已于 2026-08-06
   列出 CTNNA2、SDK1 和 TCF7L1，早于 2026-08-09 启动的本轮 v2 重算。因此，v2
   结果用于重新核验和支持三个位点，不能表述为“看完本轮新结果后才选择”。投稿
   文本采用中性且可核查的表述：三个位点被带入湿实验验证，修订分析分别报告其
   区间复发、双 caller、断点质量、TCF7L1 gene-level 稳健性和长度匹配依据；不使用
   “预先确定”，也不制造错误的因果时间顺序。
