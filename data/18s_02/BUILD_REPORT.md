# 18S-02 重建报告

日期：2026-09-23。此版本替换与 18S-01 高度重叠的旧 18S-02。

## 新数据概况

| 项目 | 结果 |
|---|---:|
| 序列条数 | 10,000 |
| 数据集内部唯一序列数 | 10,000 |
| Class 数 | 10 |
| 已分配 Family 数 | 251 |
| Genus 数 | 500 |
| 平均长度 | 378.9607 bp |
| 中位长度 | 327 bp |
| 长度范围 | 87–1,791 bp |
| 长度总体标准差 | 161.1312 bp |
| 与固定 18S-01 共享序列数 | 3,999（39.99%） |
| 相对于 18S-01 的非共享序列数 | 6,001 |
| 与旧 18S-02 共享序列数 | 4,023 |

Family 的原始标签数为 252，其中一个标签为 `NA`；792 条记录没有已分配的 Family。数据摘要的 `taxonomy_label_counts` 延续原脚本包含 `NA` 的计数方式，本文表格与 `validation.json` 的 `assigned_taxa` 排除 `NA`。

## 构建条件与方法

- 保留现有 18S-01 及其余四组数据不变。
- 保留 10,000 条序列与 10 个 Class；Genus 数经用户确认由 300 调整为 500。
- 与 18S-01 的精确序列重叠严格少于 4,000 条，即上限为 3,999。
- 使用原 `threshold_rearrangement_experiment/18S/` 的 100 对 FASTA/taxonomy 文件，共 1,000,000 条源记录；沿用原脚本的去重、分类缺失过滤和 18S Family 规范化方式，得到 57,386 条合格唯一序列。
- 同时对分类选择和属内配额分配施加共享序列硬预算；各属优先抽取不在 18S-01 中的序列。
- seed=1802；同一个 Genus 名称只选择一条分类路径。没有根据聚类效果挑选数据。

原脚本中的 `avoid_records` 只是软偏好，不能保证重复上限。新的 `--reference-fasta` 与 `--max-overlap` 按序列内容限制重叠，不依赖重新生成的序列 ID。

可行性求解确认，在当前源池与标签约束下，10,000 条、10 个 Class、300 个 Genus、共享上限 3,999 的组合不可行；500 个 Genus 的组合可行。新版本在硬预算内使用逐层配额分配尽量均衡属内数量，不进行反复的全局最小极差优化。实际每属 4–358 条，均值 20 条，总体标准差 42.5974。此分布不应描述成均匀或已达到最优平衡。

## Class 组成

| Class | 序列数 |
|---|---:|
| Aconoidasida | 774 |
| Arachnida | 1,107 |
| Bangiophyceae | 486 |
| Chlorophyceae | 1,222 |
| Conoidasida | 893 |
| Dothideomycetes | 710 |
| Eurotiomycetes | 806 |
| Insecta | 1,463 |
| Saccharomycetes | 1,408 |
| Trebouxiophyceae | 1,131 |

## 验证

已检查 FASTA/CSV 各 10,000 条、ID 一一对应、序列内容无重复、10 个 Class 和 500 个 Genus，以及重叠严格少于 4,000 条。全部输出序列、source_id 和七层分类标签均与构建时的规范化源池逐条一致。独立统计见 `validation.json`。

序列集合 SHA-256：`cffbcd08530e3b5833f6e26eda5aebae2552aecaf760af3cbb0b09a378ba5b4e`。

## 复现

在项目根目录 `/Users/pp/Documents/ChatGPT/FF` 执行：

```bash
/opt/anaconda3/bin/python -u threshold_rearrangement_experiment/build_rearrangement_datasets.py \
  --datasets 18s_02 --genera 500 \
  --reference-fasta rearrangement_clustering/data/18s_01/18s_01.fasta \
  --max-overlap 3999 \
  --output-root threshold_rearrangement_experiment/data/rebuilt_18s_02_20260923
```

新生成数据位于该输出目录的 `18s_02/` 子目录，并已复制至正式输入 `rearrangement_clustering/data/18s_02/`。原生成脚本的默认六组规格不变；重建此版本应使用以上完整参数。

## 旧结果与新版本

旧 18S-02、对应缓存及完整旧结果快照保存在 `rearrangement_clustering/archive/18s_02_before_rebuild_20260923/`。

后续实验更新已于 2026-09-23 完成：新版本的 RC 和四个基线、RP、单中心覆盖及表示统计均已重跑，六组汇总表和图已更新，其他五组输入及指标保持不变。详细结果见 `rearrangement_clustering/results/18s_02_refresh_report.md`。
