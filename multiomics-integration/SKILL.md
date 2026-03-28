---
name: multiomics-integration
description: 多组学联合分析技能包。提供转录组+蛋白组、蛋白组+乳酸化、蛋白组+磷酸化、蛋白组+代谢组、三组学修饰串扰（蛋白组+磷酸化+乳酸化）五大联合分析框架的完整分析流程、代码模板和技术指导。用于：(1) 开展多组学联合分析服务 (2) 搭建自动化分析流程 (3) 课题设计与可行性评估 (4) 生信分析方法咨询。
---

# Multi-Omics Integration Skill

多组学联合分析技能包，专注蛋白组及其修饰组学与其他组学的整合分析。

## 快速开始

### 支持的分析框架

| 框架 | 组学组合 | 核心价值 |
|------|---------|---------|
| 1 | 转录组 + 蛋白组 | 转录后调控机制 |
| 2 | 蛋白组 + 乳酸化 | 新型PTM机制研究 |
| 3 | 蛋白组 + 磷酸化 | 信号转导网络 |
| 4 | 蛋白组 + 代谢组 | 功能表型关联 |
| 5 | 蛋白组 + 磷酸化 + 乳酸化 | 修饰串扰机制 |

### 基础分析流程（所有框架通用）

```
数据输入 → 差异分析 → 重叠分析 → 关联分析 → 功能富集 → 网络可视化 → 临床关联
```

## 分析框架详情

### 框架1：转录组 + 蛋白组
**关键分析点：**
- DEGs ∩ DEPs 差异重叠（4象限分类）
- mRNA-蛋白相关性（Spearman/Pearson）
- 翻译效率 = 蛋白/mRNA 比值
- PPI网络 + 转录因子-Target网络

**参考文档：** [references/framework1_transcriptome_proteome.md](references/framework1_transcriptome_proteome.md)

### 框架2：蛋白组 + 乳酸化
**关键分析点：**
- 乳酸化位点鉴定（Kla antibody enrichment）
- 7种表达-修饰关系分类（Type A-G）
- 乳酸化motif分析 + Writer/Eraser预测
- 代谢通路（糖酵解/线粒体）富集

**参考文档：** [references/framework2_proteome_lactylation.md](references/framework2_proteome_lactylation.md)

### 框架3：蛋白组 + 磷酸化
**关键分析点：**
- 磷酸化位点定量（TiO2/IMAC enrichment）
- KSEA激酶活性推断
- 磷酸化-only蛋白（信号转导核心）
- 信号通路网络（PI3K/AKT、MAPK等）

**参考文档：** [references/framework3_proteome_phosphorylation.md](references/framework3_proteome_phosphorylation.md)

### 框架4：蛋白组 + 代谢组
**关键分析点：**
- Spearman相关性热图
- CCA典型相关分析
- 酶-底物因果调控网络
- 多组学Biomarker（LASSO/随机森林）

**参考文档：** [references/framework4_proteome_metabolome.md](references/framework4_proteome_metabolome.md)

### 框架5：三组学修饰串扰（进阶）
**关键分析点：**
- 磷酸化-乳酸化位点级串扰分析
- 双修饰蛋白分类（协同/竞争/独立）
- 激酶-乳酸化酶联合调控网络
- 通路共调控模式（糖酵解/mTORC1/免疫）

**参考文档：** [references/framework5_triple_modification.md](references/framework5_triple_modification_crosstalk.md)

## 分析代码

### Python环境依赖
```
pandas numpy scipy scikit-learn matplotlib seaborn plotly
networkx python-louvain
statsmodels
```

### 代码脚本位置
```
scripts/
├── 1_transcriptome_proteome_integration.py
├── 2_proteome_lactylation_integration.py
├── 3_proteome_phosphorylation_integration.py
├── 4_proteome_metabolome_integration.py
├── 5_triple_modification_crosstalk.py
└── common_functions.py
```

### 快速运行
```bash
# 单框架分析
python scripts/1_transcriptome_proteome_integration.py --input ./data --output ./results

# 全部框架
python run_all_frameworks.py --input ./data --output ./results
```

## 数据格式要求

### 标准输入格式
- 表达矩阵：TSV/CSV，列=样本，行=分子
- 差异分析结果：含log2FC、P-value、Regulation列
- 修饰位点数据：含位点信息（蛋白、氨基酸、位点）

### 推荐目录结构
```
project/
├── input/
│   ├── proteome.tsv
│   ├── phosphoproteome.tsv
│   ├── lactylome.tsv
│   ├── transcriptome.tsv
│   └── metabolome.tsv
├── scripts/
├── output/
└── config.yaml
```

## 参数配置

编辑 `config/pipeline_config.yaml` 调整：
- 差异分析阈值（默认 |log2FC|>1, P<0.05）
- 富集分析阈值（默认 P<0.05, Q<0.05）
- 网络构建参数（边过滤、布局算法）

## 可视化输出

| 图表类型 | 适用框架 | 说明 |
|---------|---------|------|
| Venn图 | 全部 | 差异分子重叠 |
| 热图 | 全部 | 表达/修饰相关性 |
| 火山图 | 全部 | 差异分子分布 |
| 网络图 | 1/3/5 | PPI/通路/串扰网络 |
| 富集条形图 | 全部 | GO/KEGG富集结果 |
| 串扰矩阵 | 5 | 磷酸化-乳酸化互作 |

## 临床转化应用

- **生物标志物**：多组学联合panel > 单组学
- **预后模型**：机器学习（LASSO、XGBoost、随机森林）
- **药靶发现**：激酶-底物网络核心节点

## 常见问题

**Q: 如何处理缺失值？**
A: 蛋白组数据建议用KNN或随机森林插补；代谢组用missForest

**Q: 磷酸化数据如何定量？**
A: 推荐MaxQuant + PhosphoRS；也可用pFind、Pirat

**Q: 多组学如何选择Biomarker？**
A: 先用差异交集+相关性筛选候选，再用机器学习模型验证

## 技术支持

- 分析问题：参考各框架reference文档
- 代码bug：检查依赖版本（Python>=3.8）
- 数据格式：参考示例数据（generate_sample_data.py）
