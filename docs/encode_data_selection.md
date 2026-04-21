# ENCODE Data Selection

本文件记录当前项目从 ENCODE 选择 eCLIP 数据的规则、筛选边界和可追溯依据。它用于解释为什么当前分析使用这些 experiment、file 和 primary peak sets。

## 1. 数据来源

当前数据来源为 ENCODE Project。

当前主分析对象：

- assay：`eCLIP`
- organism：`Homo sapiens`
- assembly：`GRCh38`
- file type：peak BED files

当前 seed RBP：

- `HNRNPK`
- `TARDBP`

当前 biosample：

- `HepG2`
- `K562`

当前阶段只围绕 peak 文件展开。FASTQ、BAM 和 bigWig 暂不进入主分析流程。

## 2. Experiment-level 选择规则

ENCODE experiment metadata 收集脚本：

```text
scripts/01_metadata/01_collect_encode_eclip_metadata.py
```

默认 experiment-level 过滤条件：

- `type = Experiment`
- `assay_title = eCLIP`
- `status = released`
- `organism = Homo sapiens`
- `target.label` 来自 seed RBP 列表
- `biosample_ontology.term_name` 来自 biosample 列表或命令行参数

输入控制表：

```text
metadata/rbp/rbp_seed_list.tsv
metadata/biosample/biosample_list.tsv
```

主要输出：

```text
metadata/encode/experiments.tsv
metadata/encode/files.tsv
metadata/encode/audit.tsv
```

## 3. Experiment QC 标记

`metadata/encode/experiments.tsv` 中使用以下字段记录 QC 初筛：

- `audit_warning`
- `qc_keep`
- `qc_reason`

当前规则是保守工程标记：

- released 且无明显严重 audit 的 experiment 标记为 `qc_keep = yes`。
- 出现 warning、internal action 或其他 audit flag 的 experiment 标记为 `manual_check` 或在 `qc_reason` 中记录原因。

注意：

- `qc_keep = yes` 不代表最终生物学质量认证。
- `manual_check` 不代表一定排除，而是提示后续解释时需要更谨慎。
- ENCODE audit 信息应始终保留在 `metadata/encode/audit.tsv` 中，便于回溯。

## 4. File-level 选择规则

当前主流程只使用满足以下条件的 ENCODE file：

- `output_type = peaks`
- `file_format = bed`
- `assembly = GRCh38`
- 有可用 `download_url`

这些规则的目标是让所有后续 peak 文件具有统一坐标体系和可比较格式。

暂不使用的文件类型：

- FASTQ：原始 reads，距离当前 peak-level 问题太远。
- BAM：比对结果，当前阶段不重新进行 peak calling。
- bigWig：信号轨迹，可用于后续可视化或 signal-level 分析，但不进入当前主流程。

## 5. 下载 manifest

下载 manifest 构建脚本：

```text
scripts/01_metadata/02_build_download_manifest.py
```

下载脚本：

```text
scripts/02_download/download_encode_data.py
```

下载控制表：

```text
metadata/manifest/download_manifest.tsv
```

下载状态记录：

```text
logs/download_status.tsv
logs/download_encode_data.log
```

原始 peak 文件保存路径：

```text
raw/encode/peaks/<target>/<biosample>/<experiment>/<file>.peaks.GRCh38.bed
```

约定：

- `raw/encode/peaks/` 中的文件视为原始输入，不在原地修改。
- 若下载失败，应记录在下载状态表和日志中。
- 后续处理只依赖 manifest 和本地文件状态，不直接手工挑选 raw 文件。

## 6. Peak 文件体检

peak inventory 脚本：

```text
scripts/03_processing/01_inventory_peak_files.py
```

主要输出：

```text
processed/peaks/inventory/peak_file_inventory.tsv
processed/peaks/inventory/experiment_peak_summary.tsv
logs/peak_inventory.log
```

体检重点：

- 本地文件是否存在。
- 文件大小是否合理。
- 是否能读取为 BED-like 文本。
- 数据行数和 peak 数是否合理。
- 前三列是否可作为 `chrom/start/end`。
- 坐标是否满足 `start >= 0` 和 `end > start`。

inventory 结果是 primary peak file 选择的基础。

## 7. Primary peak file 选择

主 peak 文件选择脚本：

```text
scripts/03_processing/02_select_primary_peak_files.py
```

主要输出：

```text
processed/peaks/selected/primary_peak_files.tsv
processed/peaks/selected/peak_file_ranking.tsv
logs/select_primary_peak_files.log
```

当前每个 experiment 选择一个 primary peak file。选择依据包括：

- inventory 状态是否正常。
- peak 数是否达到合理阈值。
- `output_type` 是否为 `peaks`。
- 文件大小是否合理。
- `file_format` 是否为 `bed`。
- `assembly` 是否为 `GRCh38`。

`selection_score` 用于在同一个 experiment 的多个候选 peak 文件中排序。`selection_reason` 记录每个加分项，便于回溯。

重要解释边界：

- primary peak file 是工程代表文件，不是所有 replicate 或 file 的完整统计模型。
- 后续做强生物学结论前，应考虑 replicate-aware analysis。

## 8. 当前 round 1 使用的 primary peak files

当前 round 1 和 round 2 的核心 primary peak files：

| Experiment | Target | Biosample | Primary file | Peak count |
|---|---|---|---|---:|
| `ENCSR828ZID` | `HNRNPK` | `HepG2` | `ENCFF754XAQ` | 159594 |
| `ENCSR953ZOA` | `HNRNPK` | `K562` | `ENCFF969JDI` | 182444 |
| `ENCSR187VEQ` | `TARDBP` | `HepG2` | `ENCFF626XAA` | 160498 |
| `ENCSR584TCR` | `TARDBP` | `K562` | `ENCFF606RXB` | 107504 |
| `ENCSR720BJU` | `TARDBP` | `K562` | `ENCFF381KTD` | 42044 |

这些文件构成当前小规模 exploratory dataset。它们适合用于流程验证和候选信号生成，但不足以支持普适性结论。

## 9. 标准化前后的质量检查

标准化脚本：

```text
scripts/03_processing/03_standardize_primary_peak_files.py
```

主要输出：

```text
processed/peaks/standardized/standardized_peak_file_summary.tsv
logs/standardize_primary_peak_files.log
```

标准化后应检查：

- `written_peak_count` 是否接近预期 peak 数。
- 是否存在大量 `dropped_bad_coordinates`。
- 是否存在大量 `dropped_end_le_start`。
- 是否存在大量 `dropped_negative_start`。
- 是否存在异常重复行。

如果一个文件在标准化阶段大量丢失 peaks，应回到 raw file、inventory 和 ENCODE metadata 检查原因。

## 10. 区域注释后的质量检查

区域注释脚本：

```text
scripts/04_analysis/01_annotate_peaks_with_gtf_regions.py
```

主要输出：

```text
processed/peaks/annotations/annotation_summary.tsv
processed/peaks/annotations/annotation_file_manifest.tsv
logs/annotate_peaks_with_gtf_regions.log
```

当前重要 QC 指标：

- `annotation_rate_non_intergenic`
- `chrom_compatible`
- `gtf_chrom_style`
- `peak_chrom_style`
- `status`

当前 round 1 数据中，所有 primary peak files 的 `annotation_rate_non_intergenic` 约为 99.9%。这说明 peak 文件与 GTF 注释坐标体系整体匹配良好。

## 11. 当前 priority follow-up 选择

Round 2 的 priority follow-up 来自：

```text
results/tables/exploratory_priority_comparisons.tsv
```

当前最优先比较对象：

```text
TARDBP HepG2 ENCSR187VEQ ENCFF626XAA
vs
TARDBP K562 ENCSR584TCR ENCFF606RXB
```

选择理由：

- 同一 RBP 跨 biosample 比较，解释路径相对清晰。
- 主导差异区域为 `intron`。
- sample-level absolute delta percent 最高。
- aggregated-level TARDBP HepG2 vs K562 也支持同方向信号。

该比较用于后续 shared/specific peak sets、region-specific subsets、gene mapping 和 functional input gene lists。

## 12. 扩展数据集时的更新规则

当新增 RBP 或 biosample 时，应按以下顺序更新：

1. 更新 `metadata/rbp/rbp_seed_list.tsv` 或扩展列表。
2. 更新 `metadata/biosample/biosample_list.tsv`。
3. 重新收集 ENCODE metadata。
4. 重新生成 download manifest。
5. 下载 peak 文件。
6. 运行 peak inventory。
7. 重新选择 primary peak files。
8. 标准化 primary peak files。
9. 准备 annotation inputs。
10. 重新运行 region annotation。
11. 重新汇总 region distribution。
12. 重新生成 priority comparisons。

新增数据不应手工绕过 metadata、manifest 和 inventory 三个控制步骤。

## 13. 当前选择策略的局限性

当前 ENCODE data selection 适合小规模 exploratory analysis，但存在以下限制：

- 只使用 ENCODE 已发布 peak files，没有重新从 FASTQ/BAM 做统一 peak calling。
- 每个 experiment 当前选择一个 primary peak file，尚未充分利用 replicate-level 信息。
- ENCODE audit 只做初步标记，尚未形成严格排除标准。
- 只分析 `GRCh38` peak files，跨 assembly 数据暂不合并。
- 当前 seed dataset 规模小，不能支持普适性 RBP binding atlas 结论。

这些限制应在 exploratory findings 和后续 manuscript-style 文档中明确说明。
