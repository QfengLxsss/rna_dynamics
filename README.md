# RNA Dynamics / RBP Binding Site Collection Project

## 📊 Current Results and Findings

The initial exploratory analysis (round 1) has been completed based on a small ENCODE eCLIP dataset (HNRNPK, TARDBP in HepG2 and K562).

Key observations include:
- intron-dominant binding patterns across all samples
- strong sample-dependent variation for TARDBP
- CDS preference differences between HNRNPK and TARDBP in K562

For full details, see:

👉 **[Exploratory Findings Round 1](docs/exploratory_findings_round1.md)**

## 1. 项目概述

本项目旨在围绕 RNA-binding proteins（RBP，RNA 结合蛋白）的结合位点数据，构建一套从 ENCODE 数据收集、元数据整理、peak 文件筛选、原始文件下载、文件体检、主峰文件选择、标准化输出、区域注释、区域组成比较到探索性结果总结的可复现分析流程。

当前阶段的目标不是直接做大规模结论，而是：

1. 建立一条稳定、清晰、可扩展的数据工程与预处理管线  
2. 用少量示例 RBP 和 biosample 跑通全流程  
3. 在小规模数据上识别值得进一步验证的候选信号  
4. 为后续扩展到更多 RBP、更多样本和更深入分析奠定基础  

当前项目的示例对象为：

### 示例 RBP
- `HNRNPK`
- `TARDBP`

### 示例 biosample
- `HepG2`
- `K562`

### 主要数据来源
- **ENCODE**
- 重点 assay：
  - **eCLIP**

---

## 2. 项目核心问题

本项目当前试图回答以下层层推进的问题：

1. 如何系统、可复现地从 ENCODE 中收集 RBP eCLIP 结合位点数据？
2. 如何从众多 files 中筛出适合下游分析的 peak 文件？
3. 同一个 experiment 下多个 peak 文件如何选择一个代表性主峰集？
4. 如何把原始 peak 文件标准化成可直接比较与注释的统一格式？
5. 不同 RBP 和不同 biosample 的 peak 区域组成是否存在差异？
6. 在当前小规模数据集中，哪些差异最值得后续验证与扩展？

---

## 3. 为什么围绕 peak 文件展开

### 3.1 什么是 peak 文件

在 eCLIP 等富集型测序实验中，RBP 在某些 RNA 区域会出现高于背景的结合信号。  
如果把信号沿基因组/转录组位置画出来，会形成局部富集峰，因此这些显著区域被称为 **peaks**。

peak 文件本质上是：

> 经过实验和后续 peak calling / 富集检测后得到的“显著结合区域列表”。

每一条 peak 可以理解为一个候选结合区域。

### 3.2 为什么不直接围绕 FASTQ 或 BAM

对于当前阶段的目标来说：

- FASTQ 是原始 reads，距离“结合位点结果”还太远
- BAM 是比对结果，但还不能直接告诉我们哪些区域是显著结合区域
- peak 文件已经是处理后的结果层，最接近当前真正关心的“结合位点结果”

因此，本项目当前选择以 peak 文件为核心对象。

### 3.3 peak 文件在本项目中的位置

本项目可理解为三层结构：

1. 原始证据层  
   - FASTQ  
   - BAM  

2. 结果位点层  
   - **peak files**

3. 生物学解释层  
   - peak 区域注释  
   - peak overlap  
   - motif enrichment  
   - RBP × biosample binding matrix  
   - gene-level interpretation  

当前项目已从第二层推进到第三层的起点。

---

## 4. 当前分析策略

当前采用 **小规模示范集优先、流程跑通后再扩展** 的策略。

原因如下：

1. ENCODE 的 eCLIP 数据虽然质量高，但并不是天然整理好的多组织 atlas  
2. experiment 数量和 file 类型较多，必须先建立 metadata 驱动流程  
3. 小规模示范集更容易验证流程是否稳定  
4. 当前阶段更适合先发现候选信号，再决定扩展方向

当前策略概括如下：

- 先限制在少量经典 RBP
- 先限制在 `HepG2` 和 `K562`
- 只抓 `released` 的 `human eCLIP`
- 只选择：
  - `output_type = peaks`
  - `file_format = bed`
  - `assembly = GRCh38`
- 先完成：
  - metadata 收集
  - peak 下载
  - 文件体检
  - 主峰文件选择
  - peak 标准化
  - 区域注释
  - 区域组成比较
  - 探索性结果总结

---

## 5. 当前项目目录结构

```text
rna_dynamics/
├── docs/
│   ├── project_scope.md
│   ├── data_standard.md
│   ├── encode_notes.md
│   ├── rbp_background.md
│   └── exploratory_findings_round1.md
├── logs/
├── metadata/
│   ├── biosample/
│   ├── encode/
│   ├── manifest/
│   └── rbp/
├── notebooks/
├── processed/
│   └── peaks/
│       ├── inventory/
│       ├── selected/
│       ├── standardized/
│       ├── annotation_inputs/
│       └── annotations/
├── raw/
│   ├── annotations/
│   └── encode/
│       └── peaks/
├── results/
├── scripts/
│   ├── 01_metadata/
│   ├── 02_download/
│   ├── 03_processing/
│   └── 04_analysis/
└── README.md
````

### 目录分层逻辑

* `docs/`

  * 项目设计说明、分析背景、探索性结果总结

* `metadata/`

  * 控制层
  * 存放 RBP 列表、ENCODE 元数据、biosample 信息、下载清单等

* `raw/`

  * 原始下载文件与参考注释文件
  * 原则上尽量不在此处修改原始内容

* `processed/`

  * 标准化后的 peak 文件、注释输入、注释结果、中间处理结果

* `results/`

  * 最终更偏分析解释层的结果表与汇总

* `scripts/`

  * 项目脚本，按阶段拆分

* `logs/`

  * 运行日志和状态记录

---

## 6. 当前已完成的工作流程

当前已经完成的主流程如下：

```text
→ RBP种子列表
→ ENCODE eCLIP metadata 抓取
→ 构建下载 manifest
→ 下载 peak 文件
→ peak 文件完整性与格式体检
→ 每个 experiment 选择主 peak 文件
→ 主 peak 文件标准化
→ 标准化 peak 的注释输入准备
→ 基于 GTF 的 peak 区域注释
→ 区域组成汇总与比较
→ 区域分布绘图
→ round 1 探索性结果总结
```

当前项目已经从“数据整理与预处理”正式进入“探索性结果解释阶段”。

---

## 7. 已完成脚本说明

---

### 7.1 `scripts/00_init_project_structure.py`

#### 作用

初始化项目目录结构与基础占位文件。

#### 主要功能

* 创建标准目录层级
* 创建部分 TSV 占位文件
* 创建基础文档与空目录

#### 设计意义

为整个 pipeline 提供统一骨架，避免后期目录混乱。

---

### 7.2 `scripts/01_metadata/01_collect_encode_eclip_metadata.py`

#### 作用

从 ENCODE 抓取指定 RBP、指定 biosample 的 eCLIP 元数据。

#### 输入

* `metadata/rbp/rbp_seed_list.tsv`
* 可选：`metadata/biosample/biosample_list.tsv`

#### 输出

* `metadata/encode/experiments.tsv`
* `metadata/encode/files.tsv`
* `metadata/encode/audit.tsv`

#### 当前结果

已抓取：

* `6` 个 experiments
* `132` 个 files
* `18` 条 audit 记录

---

### 7.3 `scripts/01_metadata/02_build_download_manifest.py`

#### 作用

根据 ENCODE experiments/files 元数据筛选出适合下载的 peak 文件。

#### 输入

* `metadata/encode/experiments.tsv`
* `metadata/encode/files.tsv`

#### 输出

* `metadata/manifest/download_manifest.tsv`

#### 当前筛选规则

* `qc_keep = yes`
* `output_type = peaks`
* `file_format = bed`
* `assembly = GRCh38`

#### 当前结果

从 `132` 个 file rows 中筛出：

* `15` 个 peak 文件
* `5` 个 experiments
* `2` 个 targets
* `2` 个 biosamples

---

### 7.4 `scripts/02_download/download_encode_data.py`

#### 作用

根据 manifest 批量下载 peak 文件。

#### 输入

* `metadata/manifest/download_manifest.tsv`

#### 输出

* 原始 peak 文件到 `raw/encode/peaks/...`
* `logs/download_encode_data.log`
* `logs/download_status.tsv`

#### 当前结果

* 成功下载 `15/15` 文件
* 失败 `0`
* 跳过 `0`

---

### 7.5 `scripts/02_download/download_reference_gtf.py`

#### 作用

下载参考 GTF 到项目目录结构中。

#### 输入

* 默认下载 GENCODE human release 49 basic annotation

#### 输出

* `raw/annotations/gencode.v49.basic.annotation.gtf.gz`

#### 设计意义

为后续 peak 区域注释提供参考注释文件。

---

### 7.6 `scripts/03_processing/01_inventory_peak_files.py`

#### 作用

检查下载好的 peak 文件是否存在、是否可读、格式是否一致。

#### 输入

* `metadata/manifest/download_manifest.tsv`

#### 输出

* `processed/peaks/inventory/peak_file_inventory.tsv`
* `processed/peaks/inventory/experiment_peak_summary.tsv`
* `logs/peak_inventory.log`

#### 关键修正点

* 使用 `split()` 按任意空白切分
* 按 magic bytes 自动识别 gzip，而不是仅靠文件后缀

#### 当前结果

* `15/15` 文件存在
* `15/15` 文件可读
* `15/15` 文件为一致的 `10` 列 BED-like / narrowPeak 风格

---

### 7.7 `scripts/03_processing/02_select_primary_peak_files.py`

#### 作用

为每个 experiment 选出一个代表性主 peak 文件。

#### 输入

* `metadata/encode/files.tsv`
* `metadata/manifest/download_manifest.tsv`
* `processed/peaks/inventory/peak_file_inventory.tsv`

#### 输出

* `processed/peaks/selected/primary_peak_files.tsv`
* `processed/peaks/selected/peak_file_ranking.tsv`
* `logs/select_primary_peak_files.log`

#### 当前选择策略

1. `inventory_status == ok`
2. `peak_count` 更大
3. output_type 关键词
4. 文件大小

#### 当前选出的主文件

* `ENCSR828ZID` → `ENCFF754XAQ`
* `ENCSR953ZOA` → `ENCFF969JDI`
* `ENCSR187VEQ` → `ENCFF626XAA`
* `ENCSR584TCR` → `ENCFF606RXB`
* `ENCSR720BJU` → `ENCFF381KTD`

---

### 7.8 `scripts/03_processing/03_standardize_primary_peak_files.py`

#### 作用

将选出的主 peak 文件统一标准化为 BED10 / narrowPeak 风格。

#### 输入

* `processed/peaks/selected/primary_peak_files.tsv`

#### 输出

* `processed/peaks/standardized/.../*.standardized.bed`
* `processed/peaks/standardized/standardized_peak_file_summary.tsv`
* `logs/standardize_primary_peak_files.log`

#### 当前结果

* 成功标准化 `5/5` 文件
* 无缺失
* 无报错

---

### 7.9 `scripts/03_processing/04_prepare_peak_annotation_inputs.py`

#### 作用

将标准化后的 BED10 文件转换为 annotation-ready BED6，并生成注释输入 manifest。

#### 输入

* `processed/peaks/selected/primary_peak_files.tsv`
* `processed/peaks/standardized/standardized_peak_file_summary.tsv`
* 标准化后的 BED10 文件

#### 输出

* `processed/peaks/annotation_inputs/bed6/.../*.annotation_ready.bed`
* `processed/peaks/annotation_inputs/annotation_input_manifest.tsv`
* `processed/peaks/annotation_inputs/annotation_input_summary.tsv`
* `logs/prepare_peak_annotation_inputs.log`

#### 当前结果

* 成功准备 `5/5` 个 annotation-ready BED6
* 所有文件染色体命名为 `chr_prefixed`

---

### 7.10 `scripts/04_analysis/01_annotate_peaks_with_gtf_regions.py`

#### 作用

用 GTF 为 peak 赋予主区域类别注释。

#### 输入

* `processed/peaks/annotation_inputs/annotation_input_manifest.tsv`
* `processed/peaks/annotation_inputs/bed6/.../*.annotation_ready.bed`
* `raw/annotations/gencode.v49.basic.annotation.gtf.gz`

#### 输出

* `processed/peaks/annotations/per_file/.../*.peak_region_annotations.tsv`
* `processed/peaks/annotations/annotation_summary.tsv`
* `processed/peaks/annotations/annotation_file_manifest.tsv`
* `logs/annotate_peaks_with_gtf_regions.log`

#### 注释策略

采用 **优先级注释**，即每个 peak 只分配一个主区域身份，优先级为：

```text
five_prime_utr
→ three_prime_utr
→ utr
→ CDS
→ exon
→ intron
→ gene
→ intergenic
```

#### 当前结果

* `5/5` 文件注释成功
* 无染色体命名冲突
* 所有样本非 intergenic 比例约 `99.9%`

---

### 7.11 `scripts/04_analysis/02_summarize_peak_region_distributions.py`

#### 作用

汇总 peak 区域分布，输出适合比较与绘图的多种表格。

#### 输入

* `processed/peaks/annotations/annotation_summary.tsv`

#### 输出

* `sample_region_distribution_long.tsv`
* `sample_region_distribution_wide.tsv`
* `target_biosample_region_distribution.tsv`
* `target_region_distribution.tsv`
* `biosample_region_distribution.tsv`
* `region_distribution_report.tsv`
* `logs/summarize_peak_region_distributions.log`

#### 当前结果

已成功得到 sample-level 与 aggregated-level 区域组成汇总表。

---

### 7.12 `scripts/04_analysis/03_compare_peak_region_distributions.py`

#### 作用

比较不同样本或不同 RBP 间的 peak 区域组成差异。

#### 输入

* `sample_region_distribution_wide.tsv`
* `target_biosample_region_distribution.tsv`

#### 输出

* `sample_pairwise_region_differences_long.tsv`
* `sample_pairwise_region_differences_summary.tsv`
* `aggregated_pairwise_region_differences_long.tsv`
* `aggregated_pairwise_region_differences_summary.tsv`
* `logs/compare_peak_region_distributions.log`

#### 当前关键结果

* `TARDBP: HepG2 vs K562` 的 intron 差异最大
* `HNRNPK vs TARDBP in K562` 的 CDS 差异最大

---

### 7.13 `scripts/04_analysis/04_plot_peak_region_distributions.py`

#### 作用

把区域组成和比较结果画成图。

#### 输入

* `sample_region_distribution_wide.tsv`
* `target_biosample_region_distribution.tsv`
* `sample_pairwise_region_differences_summary.tsv`
* `aggregated_pairwise_region_differences_summary.tsv`

#### 输出

目录：

* `processed/peaks/annotations/region_distribution/plots/`

主要图形：

* sample-level stacked barplot
* target+biosample stacked barplot
* sample-level heatmap
* target+biosample heatmap
* sample-level top shifts
* aggregated top shifts

#### 当前结果

全部图形已成功绘制。

---

## 8. 当前核心结果总结

截至目前，当前小规模数据集中已经出现以下明确探索性结果：

### 8.1 所有样本均以 intron peaks 为主

当前 5 个样本在区域组成上均表现为：

* intron 占主导
* CDS 次之
* intergenic 极少

### 8.2 几乎所有 peaks 都落在可解释的基因相关区域

所有样本的非 intergenic 比例约为 `99.9%`，说明：

* 注释体系匹配良好
* peak 坐标与参考注释高度兼容
* 当前流程合理可靠

### 8.3 TARDBP 在 HepG2 与 K562 间的区域偏好变化最明显

主要差异区域为：

* `intron`

### 8.4 K562 中 HNRNPK 与 TARDBP 的 CDS 偏好差异最明显

主要差异区域为：

* `CDS`

### 8.5 HepG2 中 HNRNPK 与 TARDBP 的差异较小

相比 K562，HepG2 中两者整体区域组成更接近。

---

## 9. 当前项目阶段判断

当前项目已经完成：

> **ENCODE eCLIP peak 数据收集、标准化、区域注释与第一轮探索性比较**

当前项目不再停留在“怎么收集数据”，而是已经具备：

* 可直接用于下游分析的标准化 peak 文件
* 完整的区域注释结果
* 初步可解释的探索性模式
* 一组可用于后续验证的候选信号

---

## 10. 当前结果的局限性

尽管流程已跑通，当前结果仍应被视为探索性分析，原因包括：

1. 当前仅覆盖：

   * 2 个 RBP
   * 2 个 biosample
   * 5 个主 experiment

2. 当前结果仍可能受：

   * experiment 个体差异
   * 批次效应
   * 主峰文件选择策略
   * 小样本偶然性
     的影响

3. 当前区域注释为优先级注释，适合整体区域组成比较，但不保留多重注释关系

4. 当前尚未进入：

   * difference peak set extraction
   * gene-level mapping
   * motif enrichment
   * functional enrichment
   * 大规模扩展验证

因此，当前结果最适合的定位是：

> **探索性发现（exploratory findings）与候选假设生成**

---

## 11. 下一步方向

### 优先建议：先深挖，再扩量

原因是：

* 当前小数据集已经出现了清晰信号
* 现在最重要的是把信号解释清楚
* 若此时直接大规模扩展，项目容易进入“数据越多越乱”的状态

### 建议优先深挖的两组比较

1. `TARDBP: HepG2 vs K562`
2. `HNRNPK vs TARDBP in K562`

### 推荐的探索性深挖路线

1. 定义 priority comparisons
2. 提取 comparison-specific peak sets
3. 构建 shared / specific peak 集
4. 按区域类别细分（重点关注 intron / CDS）
5. 映射到 gene 层
6. 提取候选基因
7. 做轻量功能解释

### 扩展策略建议

当上述小规模探索性分析完成后，再考虑：

* 扩展更多 RBP
* 扩展更多 biosample
* 但应围绕当前已出现的候选信号进行有针对性的扩展

---

## 12. 总结当前项目状态

本项目已经成功构建了一条以 ENCODE eCLIP 为核心、围绕示例 RBP（HNRNPK、TARDBP）与示例 biosample（HepG2、K562）展开的 **RBP 结合位点数据收集、主峰筛选、标准化、区域注释、组成比较与探索性模式识别流程**，并已识别出若干值得后续进一步验证的候选信号。

---
