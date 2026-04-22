# Exploratory Findings Round 2

**Author**：wangshuo  
**Date**：2026-04-20  
**Project**：RNA Dynamics / RBP Binding Site Collection  

---

## 1. 文档目的

本文件总结项目第二轮探索性分析的工作内容、主要结果和后续方向。

第一轮分析的重点是：基于小规模 ENCODE eCLIP 数据，跑通 peak 文件收集、标准化、区域注释和区域组成比较，并识别值得继续跟进的差异信号。

第二轮分析的重点不是增加更多 RBP 或更多细胞样本，而是围绕第一轮中最优先的比较对象，继续向下解析：

```text
-> 区域比例差异
-> 重点比较对象
-> 共享 / 特异 peak 集合
-> 按区域拆分的 peak 子集
-> 候选基因
-> 可用于功能富集的基因列表
```

因此，第二轮分析的定位是：

> 将第一轮的区域层差异，转化为可继续做功能解释的候选 peak 集合和候选基因列表。

---

## 2. 总体定位

当前第二轮结果的用途是：

- 将第一轮发现的差异信号具体化
- 识别可进一步解释的候选基因集合
- 为后续 GO / pathway 富集分析和人工文献整理准备输入
- 为后续扩展更多 RBP 或细胞样本提供优先级依据

当前结果不应被解释为普适性生物学结论，原因包括：

- 数据规模仍然较小
- 第二轮主要围绕一个重点比较对象深入
- 当前结果更适合生成假设，而不是直接作泛化结论
- 基因层面的排序尚未按基因长度、转录本复杂度或背景 peak 机会进行校正

---

## 3. 第二轮重点比较对象

第二轮主要围绕第一轮识别出的重点比较对象展开。

### 3.1 单样本层面的重点比较

当前重点比较对象为：

- 左侧：`TARDBP | HepG2 | ENCSR187VEQ | ENCFF626XAA`
- 右侧：`TARDBP | K562 | ENCSR584TCR | ENCFF606RXB`

该比较在第一轮中表现为：

- 同一 RBP 跨细胞样本差异最明显
- 主导差异区域为 `intron`
- 绝对变化幅度为 `17.1998`

### 3.2 聚合层面的支持

聚合层面的 `TARDBP | HepG2` vs `TARDBP | K562` 比较也支持这一方向：

- 主导差异区域为 `intron`
- 绝对变化幅度为 `15.9648`

因此，第二轮聚焦于：

> TARDBP 在 HepG2 与 K562 之间的样本相关 peak 结构差异。

---

## 4. 已完成的分析步骤

当前第二轮已经完成以下步骤：

```text
-> 定义重点比较对象
-> 提取重点比较对应的 peak 文件对
-> 构建共享 / 特异 peak 集合
-> 按区域类别拆分 peak 子集
-> 将 peak 子集映射到基因
-> 汇总候选基因
-> 准备用于功能富集的基因列表
```

对应脚本包括：

- `05_define_priority_comparisons.py`
- `06_extract_priority_peak_sets.py`
- `07_build_peak_overlap_sets.py`
- `08_split_peak_sets_by_region_class.py`
- `09_map_peak_sets_to_genes.py`
- `10_summarize_candidate_genes.py`
- `11_prepare_functional_gene_lists.py`

这一轮分析使项目从：

> peak 区域组成比较

推进到：

> 差异 peak 集合与候选基因层面的探索性解释。

---

## 5. 共享与特异 Peak 集合结果

基于用于注释的 BED6 文件和稳定 peak ID，当前比较被拆分为共享 peak 与样本特异 peak。

这里的左右两侧含义为：

- 左侧：`TARDBP | HepG2 | ENCSR187VEQ | ENCFF626XAA`
- 右侧：`TARDBP | K562 | ENCSR584TCR | ENCFF606RXB`

### 5.1 总体结果

| 指标 | 含义 | 数值 |
|---|---|---:|
| `left_total` | 左侧总 peak 数 | 160498 |
| `right_total` | 右侧总 peak 数 | 107504 |
| `shared_left` | 左侧中与右侧重叠的 peak 数 | 42229 |
| `shared_right` | 右侧中与左侧重叠的 peak 数 | 38330 |
| `left_specific` | 左侧特异 peak 数 | 118269 |
| `right_specific` | 右侧特异 peak 数 | 69174 |
| `pair_links` | 共享 peak 配对关系数 | 44664 |
| `left_shared_fraction` | 左侧共享比例 | 0.263112 |
| `right_shared_fraction` | 右侧共享比例 | 0.356545 |
| `left_specific_fraction` | 左侧特异比例 | 0.736888 |
| `right_specific_fraction` | 右侧特异比例 | 0.643455 |

### 5.2 结果解释

TARDBP 在 HepG2 与 K562 之间确实存在一部分共享 peak 区域：

- HepG2 侧约 26.3% 的 peaks 与 K562 侧 peaks 重叠
- K562 侧约 35.7% 的 peaks 与 HepG2 侧 peaks 重叠

但两侧仍有较大比例的样本特异 peaks：

- HepG2 特异 peaks：118269
- K562 特异 peaks：69174

这说明当前观察到的 HepG2 / K562 差异，主要需要从两侧的样本特异 peak 集合中继续解析，而不只是从共享 peaks 中解释。

---

## 6. 按区域拆分的 Peak 子集

在共享 / 特异 peak 集合基础上，进一步按 `primary_region_class` 拆分。当前第二轮最值得关注的区域类别是：

- `intron`
- `CDS`

### 6.1 关键结果

| Peak 集合 | 区域 | peak 数 | 集合内占比 |
|---|---|---:|---:|
| `left_specific` | `intron` | 95130 | 0.804353 |
| `right_specific` | `intron` | 37474 | 0.541735 |
| `shared_left` | `intron` | 33854 | 0.801677 |
| `shared_right` | `intron` | 30431 | 0.793921 |
| `right_specific` | `CDS` | 20749 | 0.299954 |
| `left_specific` | `CDS` | 7279 | 0.061546 |

### 6.2 区域层观察

**观察 1：HepG2 特异 peaks 保持 intron-heavy 特征**

`left_specific | intron` 占比约 80.4%。这说明 HepG2 特异的 TARDBP peaks 仍然高度集中在 intron 区域。

**观察 2：K562 特异 peaks 中 intron 占比降低**

`right_specific | intron` 占比约 54.2%，明显低于 HepG2 特异 peaks。这提示 K562 特异 peaks 不再像 HepG2 特异 peaks 那样主要由 intron 区域构成。

**观察 3：K562 特异 peaks 中 CDS-associated 子集比例明显升高**

`right_specific | CDS` 占比约 30.0%，而 `left_specific | CDS` 约 6.2%。这支持第一轮中观察到的 K562 侧 CDS 偏好变化，并将该变化落实到更具体的 K562 特异 CDS-associated peak 子集。

这里的 `CDS-associated` 指当前优先级注释下 `primary_region_class = CDS` 的 peak 子集，不代表这些 peaks 只与 CDS 区域发生唯一重叠。

---

## 7. 基因层面映射

在按区域拆分的 peak 子集基础上，第二轮进一步将重点 peak 集合映射到基因层面。

本轮重点关注 4 个集合：

- `left_specific | intron`
- `right_specific | intron`
- `left_specific | CDS`
- `right_specific | CDS`

### 7.1 基因映射结果

| Peak 集合 | 区域 | peak 行数 | 基因数 | 排名第一基因 | 对应 peak 命中数 |
|---|---|---:|---:|---|---:|
| `left_specific` | `intron` | 95130 | 14940 | `PTPRN2` | 1259 |
| `right_specific` | `intron` | 37474 | 8935 | `SMYD3` | 208 |
| `left_specific` | `CDS` | 7279 | 2932 | `APOB` | 119 |
| `right_specific` | `CDS` | 20749 | 5811 | `MACF1` | 56 |

### 7.2 解释边界

当前基因层面的排序主要基于 peak 命中数。它适合用于生成候选基因列表，但还没有校正：

- 基因长度
- 转录本数量
- 转录本结构复杂度
- 每个基因可被 peak 覆盖的背景机会
- 细胞样本特异的表达背景

因此，排名靠前的基因应被视为探索性候选基因，而不是已经校正后的富集基因。

---

## 8. 候选基因结果

第二轮对 4 个重点 peak 集合进行了候选基因汇总，并导出排名靠前的候选基因。

### 8.1 HepG2 特异 Intron 候选基因

`left_specific | intron`

- 排名第一基因：`PTPRN2`
- 前 5 个基因：`PTPRN2`、`LRMDA`、`PLCB1`、`SHANK2`、`MAGI1`

该集合规模最大，提示 HepG2 特异 intron peaks 对应一批覆盖面较广的 intron-associated 候选基因。

### 8.2 K562 特异 Intron 候选基因

`right_specific | intron`

- 排名第一基因：`SMYD3`
- 前 5 个基因：`SMYD3`、`PRKCB`、`MGC4859`、`CCDC26`、`ARL15`

该集合提示 K562 中也存在一批 intron-associated TARDBP 候选基因，但其前列基因组成与 HepG2 特异 intron 集合不同。

### 8.3 HepG2 特异 CDS 候选基因

`left_specific | CDS`

- 排名第一基因：`APOB`
- 前 5 个基因：`APOB`、`FN1`、`LAMA5`、`TF`、`ABCC2`

这些基因带有 HepG2 / 肝细胞背景相关线索，说明 HepG2 特异 peaks 中也存在 CDS-associated 候选基因子集，只是其比例低于 K562 特异 CDS 子集。

### 8.4 K562 特异 CDS 候选基因

`right_specific | CDS`

- 排名第一基因：`MACF1`
- 前 5 个基因：`MACF1`、`KMT2C`、`RNF213`、`UBR4`、`TPR`

该集合与 K562 特异 CDS-associated peak 子集直接对应，是第二轮中最值得优先进入功能解释和文献整理的候选集合。

---

## 9. 功能分析输入准备结果

第二轮的最后一步，是将候选基因整理为适合后续功能解释和富集分析的输入文件。

已成功导出：

- 每个重点 peak 集合的 `top_genes.tsv`
- 每个重点 peak 集合的 `gene_symbols.txt`
- 每个重点 peak 集合的 `gene_ids.txt`
- 所有选中 peak 集合的合并基因列表
- 反复出现的基因列表

当前功能分析输入文件位于：

```text
results/tables/functional_inputs/
```

当前导出的重点基因集合：

| Peak 集合 | 区域 | 导出基因数 | 排名第一基因 | 前 5 个基因 |
|---|---|---:|---|---|
| `left_specific` | `CDS` | 20 | `APOB` | `APOB;FN1;LAMA5;TF;ABCC2` |
| `left_specific` | `intron` | 20 | `PTPRN2` | `PTPRN2;LRMDA;PLCB1;SHANK2;MAGI1` |
| `right_specific` | `CDS` | 20 | `MACF1` | `MACF1;KMT2C;RNF213;UBR4;TPR` |
| `right_specific` | `intron` | 20 | `SMYD3` | `SMYD3;PRKCB;MGC4859;CCDC26;ARL15` |

说明：

- 当前每组实际导出了 20 个基因。
- 脚本支持更高的导出数量，但当前候选基因输入表每组提供的是 top 20 级别结果。
- 功能分析输入文件只是富集分析的输入，不等同于富集分析结果。

---

## 10. 第二轮核心结论

### 10.1 TARDBP 在 HepG2 / K562 间的差异主要体现在样本特异 peak 集合

共享 peaks 确实存在，但当前差异更主要地体现在大量 HepG2 特异 peaks 与 K562 特异 peaks 中。

### 10.2 HepG2 特异 peaks 更偏 intron-heavy

HepG2 特异 TARDBP peaks 中，`intron` 占比约 80.4%，并对应较大的 intron-associated 候选基因集合。

### 10.3 K562 特异 peaks 同时包含 intron 与更突出的 CDS-associated 子集

K562 特异 peaks 中，`CDS` 占比约 30.0%，明显高于 HepG2 特异 CDS 子集。这将第一轮中的 CDS 偏好变化，进一步落实到了样本特异 peak 集合与候选基因层面。

### 10.4 基因层面结果提供了下一步解释入口

当前排名靠前的候选基因包括：

- HepG2 特异 intron：`PTPRN2`
- K562 特异 intron：`SMYD3`
- HepG2 特异 CDS：`APOB`
- K562 特异 CDS：`MACF1`

这些基因是后续功能富集和文献整理的候选入口，而不是最终机制结论。

---

## 11. 科学意义

第二轮的主要价值在于，它把第一轮的区域组成差异进一步具体化。

第一轮主要回答：

- 哪些 RBP / 细胞样本比较存在较强区域比例差异？
- 哪些区域类别贡献了主要变化？

第二轮进一步回答：

- 哪些共享 / 特异 peak 集合与这些差异相关？
- 哪些按区域拆分的 peak 子集最值得继续解释？
- 这些 peak 子集对应哪些候选基因？
- 哪些基因列表可以进入下一轮功能解释？

因此，第二轮不是最终解释层，而是从区域层信号到基因层跟进分析的桥接层。

---

## 12. 当前局限性

### 12.1 数据规模有限

当前第二轮主要围绕一个重点比较对象展开：

- `TARDBP | HepG2`
- `TARDBP | K562`

因此，结果尚不能推广到所有 RBP 或所有细胞样本。

### 12.2 共享 / 特异 peak 的定义仍是简化策略

当前共享 peaks 由区间重叠定义，默认至少 1 bp 重叠。该策略适合探索性跟进，但不同重叠阈值可能影响共享 / 特异 peak 的划分。

### 12.3 基因层面排序尚未校正背景

当前排名靠前的基因主要按 peak 命中数排序。该排序可能受到基因长度、转录本数量、注释复杂度和可被 peak 覆盖的区域大小影响。因此，排名靠前的基因只能作为候选解释入口。

### 12.4 功能富集尚未完成

第二轮已经准备好可用于功能富集的基因列表，但尚未正式完成：

- GO 富集
- pathway 富集
- 背景基因集合定义
- 多重检验校正
- 系统性文献整合

这些内容应进入第三轮分析。

---

## 13. 第三轮计划方向

对第二轮产出的基因列表做轻量功能解释

优先级如下：

### 优先级 A：`right_specific | CDS`

这是当前最能体现 K562 特异 CDS-associated 信号的集合，可以优先做 GO / pathway 富集和前列基因文献整理。

### 优先级 B：`right_specific | intron`

用于理解 K562 特异 intron-associated 候选基因，并与 `right_specific | CDS` 对照。

### 优先级 C：`left_specific | intron`

用于解释 HepG2 特异 intron-heavy peak 集合，并作为 K562 特异集合的重要对照。

### 优先级 D：`left_specific | CDS`

用于提供 HepG2 特异 CDS-associated 对照集合。

第三轮计划完成：

1. 为 4 组基因列表定义合适的背景基因集合
2. 运行 GO / pathway 富集分析
3. 整理排名靠前的富集条目和候选基因
4. 对 `right_specific | CDS` 做重点文献解释

---

## 14. 总结

第二轮已经完成从区域层差异到基因层跟进分析的探索性推进。

本轮最重要的成果包括：

1. 明确 `TARDBP | HepG2 vs K562` 是当前最值得深挖的重点比较对象
2. 构建共享 / 特异 peak 集合
3. 发现 HepG2 特异 peaks 更 intron-heavy，而 K562 特异 peaks 中存在更突出的 CDS-associated 子集
4. 将重点 peak 子集映射到候选基因
5. 导出 4 组可用于功能富集的基因列表

总结：

> 第二轮将第一轮的区域组成差异进一步解析为样本特异 peak 集合、按区域拆分的候选基因和可进入第三轮功能解释的基因列表；其中 `right_specific | CDS` 是当前最值得优先解释的核心候选集合。

