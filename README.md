# RNA Dynamics / RBP Binding Site Collection

本项目围绕RBP的 ENCODE eCLIP peak 数据，建立一套从数据收集、peak 文件筛选、标准化、区域注释、区域组成比较，到差异 peak 集合拆分、候选基因汇总和功能分析输入准备的探索性分析流程。

当前项目的定位是：

> 用一个小规模、可追溯的数据集跑通 RBP eCLIP peak 分析的核心流程，并从中提出后续值得验证的候选信号。

## 📊当前状态

当前示例数据集：

- RBP：`HNRNPK`、`TARDBP`
- 细胞样本：`HepG2`、`K562`
- 实验类型：ENCODE `eCLIP`
- 物种：人类，`Homo sapiens`
- 基因组版本：`GRCh38`
- 当前选出的主 peak 文件：5 个

当前已经完成：

- 第一轮分析：peak 收集、标准化、GTF 区域注释、区域组成比较
- 第二轮分析：重点比较对象跟进、共享/特异 peak 集合、按区域拆分 peak 子集、基因映射、候选基因汇总、功能富集输入文件准备

详细结果见👉：

- [第一轮结果](docs/exploratory_findings_round1.md)
- [第二轮结果](docs/exploratory_findings_round2.md)

### 亟待解决的问题

- [ ]  **Peak 过滤：** 对peak进行过滤，去掉一些低置信的peak，目前某些基因上的Peak数异常偏高
- [ ]  **Gene 表达过滤：** 引入 RNA-seq 表达数据，因为观察发现的区域比例差异现象可能部分来自表达谱差异，而不完全是 RBP 结合偏好差异
- [ ]  **RBP 功能分类：** 后续引入更多RBP时，对RBP做分类然后根据不同的功能分类研究
- [ ]  **多 RBP 重叠 peak：** 在同一细胞系或不同细胞系中，找多个 RBP 共同结合的区域
- [ ]  **脚本工具函数重复问题：** 目前代码过于臃肿，将常用函数抽成 scripts/utils/ 或一个小包
- [ ]  **前期数据处理优化：** 研究ENCORI数据库的应用，节约前期数据处理工作花费
- [ ]  **后期方向展望：** 偏组学分析方向还是拓展到AI方向
- [ ]  **缺乏创新性：** 更像是实践复现流程或普遍性流程编写，似乎难以找到具体的创新方向


## 项目想回答的问题

当前项目围绕以下问题推进：

1. 如何系统、可复现地从 ENCODE 收集 RBP eCLIP peak 数据？
2. 如何从多个 ENCODE 文件中选择适合下游分析的主 peak 文件？
3. 如何把不同来源的 BED-like peak 文件标准化成统一格式？
4. 不同 RBP 和不同细胞样本的 peak 区域组成是否存在差异？
5. 第一轮分析识别出的重点比较对象，能否进一步拆解为共享 peak 和样本特异 peak？
6. 这些按区域拆分的 peak 子集，能否对应到可解释的候选基因？
7. 哪些基因列表最适合进入下一轮功能富集和文献解释？


## 项目目录

```text
rna_dynamics/
├── docs/
│   ├── data_standard.md
│   ├── encode_data_selection.md
│   ├── exploratory_findings_round1.md
│   └── exploratory_findings_round2.md
├── logs/
├── metadata/
│   ├── biosample/
│   ├── encode/
│   ├── manifest/
│   └── rbp/
├── notebooks/
├── processed/
│   ├── features/
│   ├── matrices/
│   └── peaks/
├── raw/
│   ├── annotations/
│   └── encode/
├── results/
│   ├── comparisons/
│   ├── enrichment/
│   ├── figures/
│   └── tables/
├── scripts/
│   ├── 01_metadata/
│   ├── 02_download/
│   ├── 03_processing/
│   └── 04_analysis/
└── README.md
```

目录含义：

- `metadata/`：控制层，保存 RBP、细胞样本、ENCODE 实验/文件信息和下载清单。
- `raw/`：原始输入层，保存下载后的 ENCODE peak 文件和 GTF 注释文件。
- `processed/`：中间处理层，保存标准化 peak、注释输入 BED、区域注释和比较用 peak 集合。
- `results/`：结果层，保存候选基因、功能分析输入和后续解释表。
- `logs/`：运行日志和状态记录。
- `docs/`：数据标准、ENCODE 数据选择规则和每一轮探索性结果。

## 文档入口

核心文档：

- [数据标准](docs/data_standard.md)：说明数据层级、坐标约定、核心表结构、区域注释、基因映射和功能输入标准。
- [ENCODE 数据选择规则](docs/encode_data_selection.md)：说明 ENCODE 实验和文件筛选规则、质控标记、主 peak 文件选择策略。
- [第一轮结果](docs/exploratory_findings_round1.md)：总结小规模 eCLIP 数据的区域组成比较和第一轮候选信号。
- [第二轮结果](docs/exploratory_findings_round2.md)：总结 `TARDBP HepG2 vs K562` 的重点跟进、peak 集合拆分和候选基因列表。

## 当前数据集

当前主 peak 文件：

| 实验编号 | RBP | 细胞样本 | 主 peak 文件 | peak 数量 |
|---|---|---|---|---:|
| `ENCSR828ZID` | `HNRNPK` | `HepG2` | `ENCFF754XAQ` | 159594 |
| `ENCSR953ZOA` | `HNRNPK` | `K562` | `ENCFF969JDI` | 182444 |
| `ENCSR187VEQ` | `TARDBP` | `HepG2` | `ENCFF626XAA` | 160498 |
| `ENCSR584TCR` | `TARDBP` | `K562` | `ENCFF606RXB` | 107504 |
| `ENCSR720BJU` | `TARDBP` | `K562` | `ENCFF381KTD` | 42044 |

这些文件构成当前第一轮和第二轮分析的小规模探索性数据集。

## 分析流程

### 元数据收集与下载

| 步骤 | 脚本 | 作用 |
|---|---|---|
| 01 | `scripts/01_metadata/01_collect_encode_eclip_metadata.py` | 收集 ENCODE eCLIP 实验和文件信息 |
| 02 | `scripts/01_metadata/02_build_download_manifest.py` | 构建 peak 文件下载清单 |
| 03 | `scripts/02_download/download_encode_data.py` | 下载 ENCODE peak 文件 |
| 04 | `scripts/02_download/download_reference_gtf.py` | 下载 GENCODE GTF 参考注释 |

### Peak 文件处理

| 步骤 | 脚本 | 作用 |
|---|---|---|
| 05 | `scripts/03_processing/01_inventory_peak_files.py` | 检查 raw peak 文件是否存在、可读、格式是否合理 |
| 06 | `scripts/03_processing/02_select_primary_peak_files.py` | 为每个 experiment 选择一个主 peak 文件 |
| 07 | `scripts/03_processing/03_standardize_primary_peak_files.py` | 将主 peak 文件标准化为 BED10-like 格式 |
| 08 | `scripts/03_processing/04_prepare_peak_annotation_inputs.py` | 转换为用于区域注释的 BED6 文件 |

### 第一轮分析

| 步骤 | 脚本 | 作用 |
|---|---|---|
| 09 | `scripts/04_analysis/01_annotate_peaks_with_gtf_regions.py` | 基于 GTF 为 peak 做优先级区域注释 |
| 10 | `scripts/04_analysis/02_summarize_peak_region_distributions.py` | 汇总 peak 区域组成 |
| 11 | `scripts/04_analysis/03_compare_peak_region_distributions.py` | 比较不同样本或不同 RBP 的区域组成差异 |
| 12 | `scripts/04_analysis/04_plot_peak_region_distributions.py` | 绘制区域组成和区域变化图 |

### 第二轮分析

| 步骤 | 脚本 | 作用 |
|---|---|---|
| 13 | `scripts/04_analysis/05_define_priority_comparisons.py` | 从区域组成差异中定义重点比较对象 |
| 14 | `scripts/04_analysis/06_extract_priority_peak_sets.py` | 提取重点比较对应的左右 peak 文件对 |
| 15 | `scripts/04_analysis/07_build_peak_overlap_sets.py` | 构建共享 peak 和样本特异 peak 集合 |
| 16 | `scripts/04_analysis/08_split_peak_sets_by_region_class.py` | 按区域类别拆分 peak 子集 |
| 17 | `scripts/04_analysis/09_map_peak_sets_to_genes.py` | 将 peak 子集映射到基因 |
| 18 | `scripts/04_analysis/10_summarize_candidate_genes.py` | 汇总候选基因 |
| 19 | `scripts/04_analysis/11_prepare_functional_gene_lists.py` | 准备用于功能富集的基因列表 |


## 关键输出文件

### 元数据与主文件选择

- `metadata/encode/experiments.tsv`
- `metadata/encode/files.tsv`
- `metadata/encode/audit.tsv`
- `metadata/manifest/download_manifest.tsv`
- `processed/peaks/inventory/peak_file_inventory.tsv`
- `processed/peaks/selected/primary_peak_files.tsv`

### 标准化 peak 与区域注释

- `processed/peaks/standardized/standardized_peak_file_summary.tsv`
- `processed/peaks/annotation_inputs/annotation_input_manifest.tsv`
- `processed/peaks/annotations/annotation_summary.tsv`
- `processed/peaks/annotations/annotation_file_manifest.tsv`

### 第一轮结果

- `processed/peaks/annotations/region_distribution/sample_region_distribution_wide.tsv`
- `processed/peaks/annotations/region_distribution/target_biosample_region_distribution.tsv`
- `processed/peaks/annotations/region_distribution/comparisons/sample_pairwise_region_differences_summary.tsv`
- `processed/peaks/annotations/region_distribution/comparisons/aggregated_pairwise_region_differences_summary.tsv`
- `processed/peaks/annotations/region_distribution/plots/`

### 第二轮结果

- `results/tables/exploratory_priority_comparisons.tsv`
- `processed/peaks/comparison_sets/comparison_peak_set_manifest.tsv`
- `processed/peaks/comparison_sets/overlaps/overlap_summary.tsv`
- `processed/peaks/comparison_sets/by_region/region_specific_peak_summary.tsv`
- `results/tables/gene_mapping/peak_set_gene_mapping_summary.tsv`
- `results/tables/candidate_genes/candidate_gene_overview.tsv`
- `results/tables/candidate_genes/top_candidate_genes_by_peak_set.tsv`
- `results/tables/functional_inputs/functional_gene_list_summary.tsv`

## 第一轮结果概览

第一轮分析完成了区域组成层面的探索。

核心观察：

- 当前 5 个主 peak 集合均以 `intron` 区域为主。
- 所有样本的非 intergenic 注释比例约为 99.9%，说明 peak 坐标与 GTF 注释体系整体匹配良好。
- `TARDBP | HepG2` vs `TARDBP | K562` 是同一 RBP 跨细胞样本差异最明显的比较，主导差异区域为 `intron`。
- `HNRNPK` vs `TARDBP` in `K562` 是同一细胞样本内不同 RBP 差异较明显的比较，主导差异区域为 `CDS`。

详见 [第一轮结果](docs/exploratory_findings_round1.md)。

## 第二轮结果概览

第二轮选择 `TARDBP | HepG2` vs `TARDBP | K562` 作为重点跟进对象。

### 重点比较对象

- 左侧：`TARDBP | HepG2 | ENCSR187VEQ | ENCFF626XAA`
- 右侧：`TARDBP | K562 | ENCSR584TCR | ENCFF606RXB`
- 主导变化区域：`intron`
- 单样本层面的绝对变化幅度：17.1998
- 聚合层面的绝对变化幅度：15.9648

### 共享与特异 Peak

| 指标 | 数值 |
|---|---:|
| `left_total` | 160498 |
| `right_total` | 107504 |
| `shared_left` | 42229 |
| `shared_right` | 38330 |
| `left_specific` | 118269 |
| `right_specific` | 69174 |
| `pair_links` | 44664 |

### 按区域拆分后的模式

- `left_specific | intron`：95130 peaks，占比 0.804353
- `right_specific | intron`：37474 peaks，占比 0.541735
- `right_specific | CDS`：20749 peaks，占比 0.299954
- `left_specific | CDS`：7279 peaks，占比 0.061546

解释：

- HepG2 特异的 TARDBP peaks 仍然明显偏向 intron。
- K562 特异的 TARDBP peaks 中出现了更突出的 CDS-associated 子集。
- 这些观察仍然属于探索性结果，需要通过更多数据和功能富集继续验证。

### 候选基因集合

| Peak 集合 | 区域 | 基因数 | Top gene | Top hits |
|---|---|---:|---|---:|
| `left_specific` | `intron` | 14940 | `PTPRN2` | 1259 |
| `right_specific` | `intron` | 8935 | `SMYD3` | 208 |
| `left_specific` | `CDS` | 2932 | `APOB` | 119 |
| `right_specific` | `CDS` | 5811 | `MACF1` | 56 |

用于功能富集的基因列表位于：

```text
results/tables/functional_inputs/
```

详见 [第二轮结果](docs/exploratory_findings_round2.md)。

## 数据标准

重要约定：

- Peak 坐标使用 BED 规则：0-based、half-open intervals。
- 当前基因组版本为 `GRCh38`。
- 当前 GTF 参考注释为 `raw/annotations/gencode.v49.basic.annotation.gtf.gz`。
- 区域注释采用优先级分配：

```text
-> five_prime_utr
-> three_prime_utr
-> utr
-> CDS
-> exon
-> intron
-> gene
-> intergenic
```

- 每个 peak 只获得一个 `primary_region_class`。
- 当前候选基因排序主要基于 peak hit count，尚未按基因长度、转录本复杂度或背景 peak 机会进行校正。

详见 [数据标准](docs/data_standard.md)。

## 当前局限性

- 数据规模仍小：2 个 RBP、2 个细胞样本、5 个主实验。
- 主 peak 文件选择是一种工程简化，还不是完整的重复样本建模。
- 区域注释是优先级注释，不保留所有重叠区域身份。
- 共享/特异 peak 集合当前基于 interval overlap 定义。
- 基因层面排序尚未按基因长度、表达量、转录本复杂度或背景机会校正。
- 功能富集分析尚未正式完成。

## 下一步方向

第三轮分析方向：

1. 为功能富集定义合适的背景基因集合。
2. 对第二轮产出的 4 组基因列表做 GO / pathway enrichment。
3. 优先解释 `right_specific | CDS`。
4. 对 `MACF1`、`SMYD3`、`PTPRN2`、`APOB` 等 top genes 做轻量文献整理。

完成第三轮后，再更有针对性地扩展更多 RBP 和更多细胞样本。

## 项目总结

本项目已经建立了一条以 metadata 为控制层、以 peak 文件为核心对象的 ENCODE eCLIP RBP 结合位点探索性分析流程。

第一轮分析跑通了数据处理、区域注释和区域组成比较。第二轮分析则沿着最强信号继续深入，拆分出样本特异 peak 集合、按区域划分的候选基因，并准备好了功能富集输入文件。

当前值得继续跟进的候选方向是：

> `TARDBP | K562 | right_specific | CDS`

这个结果应被视为假设生成线索，后续需要通过功能富集、更大数据集和独立生物学解释继续验证。
