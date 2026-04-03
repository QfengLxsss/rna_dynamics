# RNA Dynamics / RBP Binding Site Collection Project

## 1. 项目概述

本项目旨在围绕 **RNA-binding proteins** 的结合位点数据，构建一套从 **ENCODE 数据收集、元数据整理、peak 文件筛选、原始文件下载、文件体检、主峰文件选择到标准化输出** 的可复现分析流程。

当前阶段的重点是先建立一条 **稳定、清晰、可扩展的数据工程与预处理管线**，并以少量示例 RBP 为起点跑通全流程，后续再逐步扩展到更多 RBP、更多 biosample，乃至更广义的多组织体系。

本项目当前的示例对象为：

- 示例 RBP：
  - `HNRNPK`
  - `TARDBP`

- 示例 biosample：
  - `HepG2`
  - `K562`

本项目当前使用的数据源主要为：

- **ENCODE**
- 重点 assay：
  - **eCLIP**

---

## 2. 项目核心问题

本项目试图回答的根本问题是：

> 如何系统、可复现地收集和整理不同样本中的 RBP 结合位点数据，并将其转化为适合下游注释、比较与扩展分析的标准化 peak 数据集？

1. 哪些 ENCODE experiment 是相关的？
2. 哪些 files 才是我们真正想要的结合位点结果文件？
3. 同一个 experiment 下多个 peak 文件应该保留哪一个？
4. 下载下来的 peak 文件是否完整、是否可读、格式是否一致？
5. 如何把来源一致但格式细节可能不同的 peak 文件标准化为统一格式？
6. 如何保证这条流程后续能扩展到更多 RBP 和更多样本？

---

## 3. 背景原理：为什么围绕 peak 文件展开

### 3.1 什么是 peak 文件

在 eCLIP 等富集型测序实验中，RBP 在某些 RNA 区域会产生高于背景的结合信号。  
如果把信号沿基因组/转录组位置画出来，会出现局部升高的“峰（peak）”。

因此，所谓 **peak 文件**，本质上就是：

> 经过实验和后续 peak calling / 富集检测后得到的“显著结合区域列表”。

每一条 peak 通常代表：

- 一个候选结合区域
- 某个显著富集的片段
- 一个适合进一步注释和比较的区域单位

### 3.2 为什么不直接围绕 FASTQ 或 BAM

在公开数据整理项目中，原始 FASTQ 和 BAM 都非常重要，但它们不是最适合作为当前阶段的核心分析对象。

- FASTQ：原始 reads，信息最原始，但离“结合位点结果”还很远
- BAM：比对结果，能够说明 reads 落在哪，但仍不是最终的富集区域结论
- peak 文件：已经将“显著富集区域”提取出来，更接近我们真正关心的 **结合位点结果层**

因此，对于本项目这种“构建可比较的 RBP 结合位点数据集”的任务来说，peak 文件是最合适的核心对象。

### 3.3 peak 文件在本项目中的角色

本项目可以理解为三层结构：

1. 原始证据层
   - FASTQ
   - BAM

2. 结果位点层
   - **peak files**

3. 生物学解释层
   - peak 注释
   - peak overlap
   - motif enrichment
   - RBP × biosample binding matrix
   - 区域偏好分析

peak 文件正好处于中间，是将原始实验信号转化为可比较结果的桥梁。

---

## 4. 当前分析策略

本项目目前采取的是 **小规模示范集优先、流程跑通后再扩展** 的策略。

原因如下：

1. ENCODE 中 eCLIP 数据虽然质量高，但并不是天然整理好的“多组织 atlas”
2. 不同 experiment 的文件数量多、类型复杂，不能直接盲目下载
3. 必须先建立 metadata 驱动的流程，再扩展规模
4. 先从少量 RBP 和少量 biosample 入手，更容易验证流程稳定性

当前策略总结如下：

- 先确定少量经典 RBP 作为示例
- 先限制在 `HepG2` 和 `K562`
- 先抓 `released` 的 `human eCLIP` experiments
- 只选择：
  - `output_type = peaks`
  - `file_format = bed`
  - `assembly = GRCh38`
- 先完成：
  - 元数据整理
  - 文件下载
  - 文件体检
  - 主峰文件选择
  - peak 标准化

---

## 5. 当前项目目录结构

```text
rna_dynamics/
├── docs/
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
│       └── standardized/
├── raw/
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

  * 项目说明、设计原则、背景记录

* `metadata/`

  * 控制层
  * 存放 RBP 列表、ENCODE 元数据、biosample 信息、下载清单等

* `raw/`

  * 原始下载文件
  * 原则上不在这里改写原始数据

* `processed/`

  * 标准化和中间处理结果
  * 是后续分析的主要输入来源

* `results/`

  * 下游分析结果、图表、统计表格

* `scripts/`

  * 项目脚本
  * 按流程阶段拆分

* `logs/`

  * 运行日志和状态记录

---

## 6. 已完成的工作流程总览

当前已经完成的主流程如下：

```text
RBP种子列表
→ ENCODE eCLIP metadata 抓取
→ 构建下载 manifest
→ 下载 peak 文件
→ peak 文件完整性与格式体检
→ 每个 experiment 选择一个主 peak 文件
→ 主 peak 文件标准化
```

这一流程目前已经成功跑通。

---

## 7. 已完成脚本说明

以下为当前项目中已经实现并运行的脚本及其作用。

---

## 7.1 `scripts/00_init_project_structure.py`

### 作用

初始化项目目录结构与基础占位文件。

### 主要功能

* 创建项目的标准目录层级
* 创建元数据 TSV 占位文件
* 创建部分文档与空目录

### 设计意义

该脚本的目的不是分析数据，而是保证项目从一开始就按统一结构组织，避免后期目录混乱。

---

## 7.2 `scripts/01_metadata/01_collect_encode_eclip_metadata.py`

### 作用

从 ENCODE 抓取指定 RBP、指定 biosample 的 eCLIP 元数据。

### 输入

* `metadata/rbp/rbp_seed_list.tsv`
* 可选：`metadata/biosample/biosample_list.tsv`

### 输出

* `metadata/encode/experiments.tsv`
* `metadata/encode/files.tsv`
* `metadata/encode/audit.tsv`

### 工作逻辑

1. 读取示例 RBP 列表
2. 读取示例 biosample 列表（若提供）
3. 向 ENCODE API 发起搜索请求，限定：

   * `assay_title = eCLIP`
   * `organism = Homo sapiens`
   * `status = released`
4. 抓取 experiment 级元数据
5. 展开 files 信息
6. 抓取审计信息
7. 写出为 TSV 文件

### 当前结果

当前已成功抓取到：

* `6` 个 experiments
* `132` 个 files
* `18` 条 audit 记录

---

## 7.3 `scripts/01_metadata/02_build_download_manifest.py`

### 作用

根据 ENCODE experiment/file 元数据，筛选出真正要下载的 peak 文件。

### 输入

* `metadata/encode/experiments.tsv`
* `metadata/encode/files.tsv`

### 输出

* `metadata/manifest/download_manifest.tsv`

### 工作逻辑

按保守规则筛选：

* `qc_keep == yes`
* `output_type == peaks`
* `file_format == bed`
* `assembly == GRCh38`

同时生成本地保存路径 `local_path`。

### 当前结果

从 `132` 个 file rows 中筛出：

* `15` 个可下载 peak 文件
* 涉及 `5` 个 experiments
* 涉及 `2` 个 targets
* 涉及 `2` 个 biosamples

---

## 7.4 `scripts/02_download/download_encode_data.py`

### 作用

根据 manifest 批量下载 peak 文件。

### 输入

* `metadata/manifest/download_manifest.tsv`

### 输出

* 原始下载文件落到 `raw/encode/peaks/...`
* 日志：`logs/download_encode_data.log`
* 状态表：`logs/download_status.tsv`

### 工作逻辑

1. 读取 manifest
2. 只处理 `download_flag == yes` 的条目
3. 根据 `download_url` 下载文件
4. 保存到 `local_path`
5. 支持已存在跳过、失败重试、日志记录

### 当前结果

已成功下载：

* `15/15` 个文件
* `0` 个失败
* `0` 个跳过

---

## 7.5 `scripts/03_processing/01_inventory_peak_files.py`

### 作用

检查下载好的 peak 文件是否存在、是否可读、是否格式一致。

### 输入

* `metadata/manifest/download_manifest.tsv`

### 输出

* `processed/peaks/inventory/peak_file_inventory.tsv`
* `processed/peaks/inventory/experiment_peak_summary.tsv`
* `logs/peak_inventory.log`

### 工作逻辑

1. 检查 manifest 中列出的本地文件是否存在
2. 自动识别 gzip 文件（按 magic bytes，而不是只看后缀）
3. 读取前几条有效数据行
4. 统计列数是否一致
5. 汇总到 file-level 和 experiment-level 表

### 当前结果

最终确认：

* `15/15` 文件存在
* `15/15` 文件可读
* `15/15` 文件都是一致的 `10` 列 BED-like/narrowPeak 风格

---

## 7.6 `scripts/03_processing/02_select_primary_peak_files.py`

### 作用

在每个 experiment 对应的多个 peak 文件中，选择一个代表性的主 peak 文件。

### 输入

* `metadata/encode/files.tsv`
* `metadata/manifest/download_manifest.tsv`
* `processed/peaks/inventory/peak_file_inventory.tsv`

### 输出

* `processed/peaks/selected/primary_peak_files.tsv`
* `processed/peaks/selected/peak_file_ranking.tsv`
* `logs/select_primary_peak_files.log`

### 设计背景

每个 experiment 不只对应一个 peak 文件。
如果不做筛选，后续分析会把 replicate 文件、不同处理版本混在一起，不利于标准化比较。

### 初始策略

最开始采用：

* inventory 是否正常
* output_type 关键词
* 文件大小

但后发现当前这批数据中 `output_type` 区分度不够，都只是 `peaks`。

### 升级策略

当前版本改为：

1. `inventory_status == ok`
2. **peak_count 更大（最关键）**
3. output_type 关键词
4. 文件大小

脚本中会直接读取每个 peak 文件并统计行数（即 peak 数量）。

### 当前选择结果

为 `5` 个 experiment 各选择了一个主文件：

* `ENCSR828ZID` → `ENCFF754XAQ`
* `ENCSR953ZOA` → `ENCFF969JDI`
* `ENCSR187VEQ` → `ENCFF626XAA`
* `ENCSR584TCR` → `ENCFF606RXB`
* `ENCSR720BJU` → `ENCFF381KTD`

对应 peak 数量分别为：

* `159594`
* `182444`
* `160498`
* `107504`
* `42044`

---

## 7.7 `scripts/03_processing/03_standardize_primary_peak_files.py`

### 作用

将选出的主 peak 文件统一标准化为 narrowPeak 风格的 BED10 文件。

### 输入

* `processed/peaks/selected/primary_peak_files.tsv`

### 输出

* `processed/peaks/standardized/<target>/<biosample>/<experiment>/<file>.standardized.bed`
* `processed/peaks/standardized/standardized_peak_file_summary.tsv`
* `logs/standardize_primary_peak_files.log`

### 标准化逻辑

对于每个主文件：

1. 自动识别是否为 gzip
2. 去掉空行和注释行
3. 按任意空白字符切分字段
4. 少于 3 列则丢弃
5. 若小于 10 列则补齐到 10 列
6. 若大于 10 列则只保留前 10 列
7. 坐标转整数
8. 过滤：

   * 坐标非法
   * `end <= start`
   * `start < 0`
9. 按 `(chrom, start, end, name)` 去重
10. 按坐标排序
11. 输出为统一 BED10 文件

### 当前结果

已成功标准化：

* `5/5` 个主文件
* 无错误
* 无缺失

标准化后各文件保留 peak 数如下：

* `ENCSR828ZID` → `159594`
* `ENCSR953ZOA` → `182444`
* `ENCSR187VEQ` → `160498`
* `ENCSR584TCR` → `107504`
* `ENCSR720BJU` → `42044`

---

## 8. 当前项目进度总结

到目前为止，本项目已经完成：

### 8.1 元数据层

* 收集了 ENCODE eCLIP experiment 元数据
* 收集了 file 元数据
* 收集了审计信息

### 8.2 文件控制层

* 构建了下载 manifest
* 明确筛出了高质量 GRCh38 peaks 文件

### 8.3 原始数据层

* 已完成原始 peak 文件下载

### 8.4 文件体检层

* 已检查文件存在性
* 已确认文件可读
* 已确认格式一致（BED10-like）

### 8.5 主文件选择层

* 每个 experiment 已选出一个代表性的主 peak 文件

### 8.6 标准化层

* 所有主 peak 文件均已统一标准化输出

因此，当前阶段可以定义为：

> **ENCODE eCLIP peak 数据收集与标准化预处理阶段已完成**

---

## 9. 当前已有核心结果文件

### 元数据层

* `metadata/encode/experiments.tsv`
* `metadata/encode/files.tsv`
* `metadata/encode/audit.tsv`

### 下载控制层

* `metadata/manifest/download_manifest.tsv`

### 文件体检层

* `processed/peaks/inventory/peak_file_inventory.tsv`
* `processed/peaks/inventory/experiment_peak_summary.tsv`

### 主文件选择层

* `processed/peaks/selected/primary_peak_files.tsv`
* `processed/peaks/selected/peak_file_ranking.tsv`

### 标准化结果层

* `processed/peaks/standardized/standardized_peak_file_summary.tsv`
* `processed/peaks/standardized/.../*.standardized.bed`

---

## 10. 当前阶段的科学含义

当前得到的标准化 peak 文件已经可以用于：

* 区域注释（5'UTR / CDS / 3'UTR / intron / intergenic）
* 样本间 overlap 比较
* RBP 特异结合模式分析
* motif enrichment
* 后续的 binding matrix 构建

---

## 11. 当前阶段的局限性

1. 当前仅覆盖两个示例 RBP：

   * `HNRNPK`
   * `TARDBP`

2. 当前仅覆盖两个示例 biosample：

   * `HepG2`
   * `K562`

3. 当前“主 peak 文件”的选择是基于启发式规则：

   * inventory 状态
   * peak 数量
   * output_type 关键词
   * 文件大小

4. 当前还没有正式引入基因注释层（例如 GTF/Gencode）

5. 当前还没有进入：

   * peak 注释
   * overlap 比较
   * motif 分析
   * 网络与功能分析

---

## 12. 下一步计划

### 12.1 进行 peak 注释

目标：

* 将每个标准化 peak 注释到：

  * 5'UTR
  * CDS
  * 3'UTR
  * intron
  * ncRNA
  * intergenic

意义：

* 可直接解释不同 RBP 的区域偏好
* 可比较 HepG2 与 K562 中区域分布差异

### 12.2 进行样本间 overlap 比较

优先做：

* `HNRNPK: HepG2 vs K562`
* `TARDBP: HepG2 vs K562`

意义：

* 查看不同样本中共同结合区域与特异区域
* 为“多样本 binding pattern”提供第一批定量结果

### 12.3 构建扩展框架

待小规模流程验证完成后，扩展：

* 更多 RBP
* 更多 biosample
* 更多组织来源
* 可考虑整合外部数据库

---
