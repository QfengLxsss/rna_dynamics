# Data Standard

本文件定义本项目当前阶段使用的数据约定、文件层级、核心表结构和可追溯规则。它的目标不是记录每一次分析结论，而是保证从 ENCODE eCLIP peak 数据到区域注释、peak set 分解、gene mapping 和 functional input 的过程中，所有中间数据都有稳定含义。

## 1. 数据层级

项目数据按用途分为五层：

- `metadata/`：控制层，记录 RBP、biosample、ENCODE experiment/file metadata、下载清单和筛选依据。
- `raw/`：原始输入层，保存从外部来源下载的文件。原则上不在此层修改文件内容。
- `processed/`：标准化与中间结果层，保存可复用的标准 peak、annotation-ready BED、区域注释和 comparison peak sets。
- `results/`：解释层，保存面向分析结论的汇总表、候选基因表和 functional input gene lists。
- `logs/`：运行记录层，保存脚本执行日志、状态统计和排错信息。

除 `raw/` 外，其余层级的文件可以由脚本重新生成。若后续引入 workflow manager，应优先把 `processed/` 和 `results/` 的依赖关系显式化。

## 2. 坐标与基因组版本

- 当前基因组 assembly：`GRCh38`。
- 当前参考注释：`raw/annotations/gencode.v49.basic.annotation.gtf.gz`。
- peak 坐标遵循 BED 约定：0-based、half-open interval，即 `[start, end)`。
- GTF 坐标在注释脚本内部转换为 BED-like 0-based half-open intervals 后再与 peak 比较。
- peak 与 GTF 的 chromosome naming style 必须兼容。当前流程默认检查 `chr` 前缀风格，不兼容时应停止或人工检查。
- 所有 peak 必须满足 `start >= 0` 且 `end > start`。

## 3. 文件格式通用约定

- 表格文件使用 UTF-8 编码的 TSV。
- TSV 第一行为 header。
- 缺失值以空字符串表示；不要混用 `NA`、`None`、`.` 表示表格缺失。
- 文件路径字段优先保存相对于项目根目录的路径，便于迁移和复现。
- 状态字段统一使用短标签，例如 `ok`、`error`、`missing_input`、`skipped`、`manual_check`。
- 解释性字段优先使用 `note`、`reason` 或 `selection_reason`，用于记录可追溯原因。

## 4. 核心 metadata 表

### `metadata/rbp/rbp_seed_list.tsv`

当前分析种子 RBP 列表。核心字段：

- `rbp_symbol`：RBP gene symbol，例如 `HNRNPK`、`TARDBP`。

### `metadata/rbp/rbp_expanded_list.tsv`

扩展 RBP 候选列表。用于后续从 seed set 扩大到更多 RBP。

### `metadata/biosample/biosample_list.tsv`

当前纳入分析的 biosample 列表。核心字段：

- `biosample_term_name`：ENCODE biosample 名称，例如 `HepG2`、`K562`。

### `metadata/biosample/tissue_mapping.tsv`

biosample 到组织、细胞类型或分析分组的人工映射表。用于跨 biosample 聚合和解释。

### `metadata/encode/experiments.tsv`

ENCODE experiment-level metadata。核心字段包括：

- `experiment_accession`
- `target_label`
- `assay_title`
- `biosample_term_name`
- `biosample_type`
- `organism`
- `lab`
- `status`
- `audit_warning`
- `qc_keep`
- `qc_reason`

`qc_keep` 是当前项目的保守标记，不等同于最终生物学质量结论。带有 ENCODE audit warning 的 experiment 应优先进入人工检查。

### `metadata/encode/files.tsv`

ENCODE file-level metadata。核心字段包括：

- `file_accession`
- `experiment_accession`
- `file_format`
- `output_type`
- `assembly`
- `biological_replicates`
- `technical_replicates`
- `download_url`

当前主流程只使用 `output_type = peaks`、`file_format = bed`、`assembly = GRCh38` 的 peak 文件。

### `metadata/encode/audit.tsv`

ENCODE audit 信息展开表，用于追踪 warning、internal action 或其他 QC 相关提示。

### `metadata/manifest/download_manifest.tsv`

下载控制表。每一行对应一个待下载文件，并包含本地保存路径、来源 URL 和相关 experiment/file metadata。

## 5. Raw peak 文件

原始 peak 文件保存在：

```text
raw/encode/peaks/<target>/<biosample>/<experiment>/<file>.peaks.GRCh38.bed
```

约定：

- 文件内容保持 ENCODE 下载后的原始状态。
- 不在 `raw/` 层做去重、排序、坐标修正或字段补齐。
- 文件是否可用由 inventory、selection 和 standardization 步骤记录。

## 6. Primary peak file 选择

主 peak 文件选择结果保存在：

```text
processed/peaks/selected/primary_peak_files.tsv
```

核心字段：

- `experiment_accession`
- `target_label`
- `biosample_term_name`
- `selected_file_accession`
- `selected_output_type`
- `selected_file_format`
- `selected_assembly`
- `selected_local_path`
- `selected_file_size_bytes`
- `selected_peak_count`
- `selection_score`
- `selection_reason`

`selection_score` 是工程化筛选分数，用于在同一个 experiment 的多个 peak 文件中选择一个代表性主 peak 文件。它不应被解释为生物学效应强度。

## 7. Standardized BED10

标准化 peak 文件保存在：

```text
processed/peaks/standardized/<target>/<biosample>/<experiment>/<file>.standardized.bed
```

标准输出为 BED10-like 格式，字段顺序固定：

1. `chrom`
2. `start`
3. `end`
4. `name`
5. `score`
6. `strand`
7. `signalValue`
8. `pValue`
9. `qValue`
10. `peak`

标准化规则：

- 跳过空行、comment、`track` 和 `browser` 行。
- 至少需要 3 列坐标字段。
- 少于 10 列时补齐默认值。
- 多于 10 列时保留前 10 列。
- 强制检查 `start` 和 `end` 为整数。
- 丢弃 `start < 0` 或 `end <= start` 的记录。
- 以 `(chrom, start, end, name)` 去重。
- 输出按 chromosome、start、end、name 排序。

标准化汇总表：

```text
processed/peaks/standardized/standardized_peak_file_summary.tsv
```

该表记录每个文件的原始行数、写出 peak 数、丢弃原因和运行状态。

## 8. Annotation-ready BED6

注释输入文件保存在：

```text
processed/peaks/annotation_inputs/bed6/<target>/<biosample>/<experiment>/<file>.annotation_ready.bed
```

BED6 字段顺序：

1. `chrom`
2. `start`
3. `end`
4. `peak_id`
5. `score`
6. `strand`

`peak_id` 必须稳定且唯一，用于连接 BED、region annotation、overlap set 和 gene mapping。

对应 manifest：

```text
processed/peaks/annotation_inputs/annotation_input_manifest.tsv
processed/peaks/annotation_inputs/annotation_input_summary.tsv
```

## 9. Region annotation 标准

peak 区域注释结果保存在：

```text
processed/peaks/annotations/per_file/<target>/<biosample>/<experiment>/<file>.peak_region_annotations.tsv
```

核心字段：

- `peak_id`
- `chrom`
- `start`
- `end`
- `score`
- `strand`
- `target_label`
- `biosample_term_name`
- `experiment_accession`
- `selected_file_accession`
- `primary_region_class`
- `overlap_gene_ids`
- `overlap_gene_names`
- `overlap_transcript_ids`
- `n_primary_region_features`

当前采用 priority-based annotation。若一个 peak 同时重叠多个区域类别，只赋予一个 `primary_region_class`。

区域优先级顺序固定为：

```text
five_prime_utr
three_prime_utr
utr
CDS
exon
intron
gene
intergenic
```

该策略适合做总体区域组成比较，但不保留一个 peak 的全部多重区域身份。

区域注释汇总表：

```text
processed/peaks/annotations/annotation_summary.tsv
processed/peaks/annotations/annotation_file_manifest.tsv
```

`annotation_rate_non_intergenic` 是重要 QC 指标，用于检查 peak 与 GTF 坐标体系是否整体匹配。

## 10. Region distribution 结果表

区域组成结果保存在：

```text
processed/peaks/annotations/region_distribution/
```

主要表包括：

- `sample_region_distribution_long.tsv`
- `sample_region_distribution_wide.tsv`
- `target_biosample_region_distribution.tsv`
- `target_region_distribution.tsv`
- `biosample_region_distribution.tsv`
- `region_distribution_report.tsv`

pairwise 比较结果保存在：

```text
processed/peaks/annotations/region_distribution/comparisons/
```

主要表包括：

- `sample_pairwise_region_differences_long.tsv`
- `sample_pairwise_region_differences_summary.tsv`
- `aggregated_pairwise_region_differences_long.tsv`
- `aggregated_pairwise_region_differences_summary.tsv`

这些表是 round 1 exploratory findings 的主要数据来源。

## 11. Priority comparison 与 peak set 标准

Priority comparison 表：

```text
results/tables/exploratory_priority_comparisons.tsv
```

该表把 region distribution comparison 转换为后续分析优先级。核心字段包括：

- `priority_rank`
- `priority_level`
- `priority_score`
- `comparison_scope`
- `comparison_type`
- `comparison_id`
- left/right target、biosample、experiment、file 信息
- `dominant_shift_region`
- `dominant_shift_delta_percent`
- `dominant_shift_abs_delta_percent`
- `recommended_followup_focus`
- `recommended_peak_subset`

Priority peak set manifest：

```text
processed/peaks/comparison_sets/comparison_peak_set_manifest.tsv
```

该表记录 priority comparison 对应的左右 peak 文件、annotation 文件和 materialized comparison inputs。

## 12. Overlap 与 region-specific peak sets

Overlap 结果保存在：

```text
processed/peaks/comparison_sets/overlaps/
```

主要文件：

- `overlap_summary.tsv`
- `overlap_file_manifest.tsv`
- `<comparison_id>/shared_peak_pairs.tsv`
- `<comparison_id>/shared_left.annotation_ready.bed`
- `<comparison_id>/shared_right.annotation_ready.bed`
- `<comparison_id>/left_specific.annotation_ready.bed`
- `<comparison_id>/right_specific.annotation_ready.bed`

默认 shared peak 定义为至少 1 bp interval overlap。若使用 reciprocal overlap threshold，必须在日志和结果说明中记录阈值。

Region-specific peak set 结果保存在：

```text
processed/peaks/comparison_sets/by_region/
```

主要文件：

- `region_specific_peak_summary.tsv`
- `region_specific_file_manifest.tsv`
- `<comparison_id>/<set_name>.<region_class>.annotation_ready.bed`
- `<comparison_id>/<set_name>.<region_class>.annotations.tsv`

`set_name` 当前包括：

- `shared_left`
- `shared_right`
- `left_specific`
- `right_specific`

## 13. Gene mapping 与 candidate gene 表

Peak set 到 gene 的映射结果保存在：

```text
results/tables/gene_mapping/
```

主要文件：

- `peak_set_gene_mapping_manifest.tsv`
- `peak_set_gene_mapping_summary.tsv`
- `per_set/<comparison_id>/<set_name>.<region_class>.peak_to_gene.tsv`
- `per_set/<comparison_id>/<set_name>.<region_class>.gene_summary.tsv`

gene-level candidate 表保存在：

```text
results/tables/candidate_genes/
```

主要文件：

- `candidate_gene_overview.tsv`
- `top_candidate_genes_by_peak_set.tsv`
- `candidate_gene_recurrence.tsv`
- `candidate_gene_set_summary.tsv`

这些表用于从 peak-level 差异进入 gene-level exploratory interpretation。

## 14. Functional input gene lists

功能分析输入保存在：

```text
results/tables/functional_inputs/
```

主要文件：

- `functional_gene_list_manifest.tsv`
- `functional_gene_list_summary.tsv`
- `per_set/<comparison_id>/<set_name>.<region_class>.top_genes.tsv`
- `per_set/<comparison_id>/<set_name>.<region_class>.gene_symbols.txt`
- `per_set/<comparison_id>/<set_name>.<region_class>.gene_ids.txt`
- `combined/selected_sets_union.gene_symbols.txt`
- `combined/selected_sets_union.gene_ids.txt`
- `combined/recurrent_genes_minN.tsv`

约定：

- `top_genes.tsv` 保留原始 `gene_id`，包括 Ensembl version。
- `gene_ids.txt` 使用去除 version 后的 Ensembl ID，便于外部 enrichment 工具识别。
- `gene_symbols.txt` 优先使用 gene symbol；若 gene symbol 缺失，则使用去除 version 后的 Ensembl ID 作为 fallback label。
- 当前 functional input 是富集分析的输入，不等同于富集分析结果。

## 15. 命名与可追溯性

所有脚本输出应尽量保留以下追踪字段：

- `target_label`
- `biosample_term_name`
- `experiment_accession`
- `selected_file_accession`
- `comparison_id`
- `set_name`
- `region_class`
- `status`
- `note`

路径命名中不安全字符应替换为 `_`。推荐路径结构为：

```text
<target>/<biosample>/<experiment>/<file>
<comparison_id>/<set_name>.<region_class>.<suffix>
```

任何最终结果表都应能通过这些字段追溯回原始 ENCODE file accession。

## 16. 解释边界

当前数据标准服务于 exploratory analysis。以下内容不应仅凭当前表格直接作最终结论：

- RBP 的普适性结合规律
- biosample 间的稳定生物学差异
- candidate gene 的功能因果关系
- pathway enrichment 的显著性解释

这些问题需要后续通过更多 RBP、更多 biosample、replicate-aware 分析、背景基因集定义、功能富集和独立验证继续推进。
