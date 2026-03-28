# Multi-Omics Integration Skill

一套完整的生信多组学联合分析流程，支持转录组、蛋白质组、磷酸化组、乳酸化组和代谢组的整合分析。

## 框架

| 框架 | 组学组合 | 核心分析内容 |
|------|---------|-------------|
| Framework 1 | 转录组 + 蛋白质组 | 差异重叠分析、mRNA-蛋白相关性、翻译效率分析、PPI网络 |
| Framework 2 | 蛋白质组 + 乳酸化组 | 乳酸化位点鉴定、Type A-G分类、Motif分析、代谢通路富集 |
| Framework 3 | 蛋白质组 + 磷酸化组 | 磷酸化位点定量、KSEA激酶活性推断、信号通路网络 |
| Framework 4 | 蛋白质组 + 代谢组 | 相关性热图、CCA分析、因果网络、Biomarker发现 |
| Framework 5 | 三组学修饰串扰 | 磷酸化-乳酸化位点级串扰、酶-底物网络、激酶/乳酸化酶活性 |

## 安装

```bash
pip install -r requirements.txt
```

## 数据准备

在 `input/` 目录下放置数据文件（按框架分目录）：

### Framework 1: 转录组 + 蛋白质组
```
input/framework1_transcriptome_proteome/
├── DEGs.csv              # 差异基因 (GeneSymbol, log2FC, PValue)
├── DEPs.csv              # 差异蛋白 (ProteinID, GeneSymbol, log2FC, PValue)
├── TPM.csv               # 转录组表达量 (GeneSymbol × samples)
└── protein_intensity.csv # 蛋白表达量 (ProteinID × samples)
```

### Framework 2: 蛋白质组 + 乳酸化
```
input/framework2_proteome_lactylation/
├── lactylation_sites.csv     # 乳酸化位点 (ProteinID, Site, log2FC, PValue)
└── proteome_quantification.csv # 全局蛋白定量 (ProteinID × samples)
```

### Framework 3: 蛋白质组 + 磷酸化
```
input/framework3_proteome_phosphorylation/
├── phosphorylation_sites.csv  # 磷酸化位点 (ProteinID, Site, log2FC, PValue)
└── proteome_quantification.csv
```

### Framework 4: 蛋白质组 + 代谢组
```
input/framework4_proteome_metabolome/
├── proteome_quantification.csv    # 蛋白表达量
└── metabolome_quantification.csv  # 代谢物表达量 (Metabolite × samples)
```

### Framework 5: 三组学修饰串扰
```
input/framework5_triple_modification/
├── proteome_quantification.csv    # 全局蛋白定量
├── phosphorylation_sites.csv      # 磷酸化位点
└── lactylation_sites.csv         # 乳酸化位点
```

### 样本信息（通用）
```
input/sample_info.csv
Sample,Group
Ctrl_1,Control
Ctrl_2,Control
Ctrl_3,Control
Treat_1,Treatment
Treat_2,Treatment
Treat_3,Treatment
```

## 运行

```bash
# 运行单个框架
python run_pipeline.py --framework 1

# 运行所有框架
python run_pipeline.py --framework all

# 指定输入/输出目录
python run_pipeline.py --framework 2 --input ./data --output ./results
```

## 输出

```
results/
├── 1_transcriptome_proteome/
│   ├── overlap_analysis.csv      # DEGs ∩ DEPs 重叠分析
│   ├── mRNA_protein_correlation.csv # mRNA-蛋白相关性
│   ├── enrichment_network.png    # GO/KEGG富集网络
│   ├── PPI_network.csv           # PPI网络数据
│   ├── translation_efficiency.csv # 翻译效率分析
│   ├── venn.png                  # Venn图
│   └── heatmap.png               # 热图
├── 2_proteome_lactylation/
│   ├── site_analysis.csv         # 乳酸化位点分析
│   ├── type_classification.csv   # Type A-G分类
│   └── pathway_enrichment.png    # 通路富集
├── 3_proteome_phosphorylation/
│   ├── site_quantification.csv   # 磷酸化位点定量
│   ├── ksea_results.csv          # KSEA激酶活性
│   └── signaling_network.png     # 信号通路网络
├── 4_proteome_metabolome/
│   ├── correlation_matrix.csv    # 相关性矩阵
│   ├── cca_results.csv           # CCA分析结果
│   └── biomarker_candidates.csv   # Biomarker候选
└── 5_triple_modification/
    ├── crosstalk_matrix.csv      # 串扰矩阵
    ├── enzyme_network.csv        # 酶-底物网络
    └── pathway_crosstalk.png     # 通路串扰图
```

## 配置

编辑 `config/pipeline_config.yaml` 修改分析参数：

```yaml
[differential]
fc_threshold = 1.5      # 差异倍数阈值
pvalue_threshold = 0.05 # P值阈值

[correlation]
method = "spearman"     # 相关性计算方法
threshold = 0.6         # 相关性阈值

[enrichment]
pvalue_cutoff = 0.05     # 富集分析P值截断
min_overlap = 3         # 最小重叠基因数
```

## 各框架核心分析

### Framework 1: 转录组 + 蛋白质组
- **差异重叠分析** (DEGs ∩ DEPs)
- Venn图展示重叠，分为：一致上调、一致下调、不一致
- **mRNA-蛋白表达相关性**
- Pearson/Spearman相关性，强相关/弱相关/负相关分类
- **GO/KEGG联合富集**
- 通路网络可视化
- **PPI网络构建**
- STRING数据库查询，Hub蛋白识别
- **翻译效率分析**
- TE = Protein / mRNA，翻译增强/抑制分类

### Framework 2: 蛋白质组 + 乳酸化
- **乳酸化位点鉴定 + 蛋白定量关联**
- 位点水平分析和蛋白水平分析
- **乳酸化-表达关系分类 (Type A-G)**
  - Type A: 蛋白↑ + 乳酸化↑
  - Type B: 蛋白↓ + 乳酸化↓
  - Type C: 蛋白↑ + 乳酸化↓
  - Type D: 蛋白↓ + 乳酸化↑
  - Type E: 仅蛋白变化
  - Type F: 仅乳酸化变化
  - Type G: 无显著变化
- **乳酸化Motif分析**
- 序列logos，位置频率矩阵
- **功能富集**
- 代谢通路、线粒体功能

### Framework 3: 蛋白质组 + 磷酸化
- **磷酸化位点定量 + 全局蛋白关联**
- 位点-蛋白关联表，分类：Phospho-up/down
- **KSEA激酶活性推断**
- Z-score计算，激酶活性状态判定
- **信号通路网络构建**
- 激酶-底物网络，通路图绘制
- **磷酸化-only蛋白分析**
- 信号转导核心节点，磷酸化的独立效应

### Framework 4: 蛋白质组 + 代谢组
- **相关性热图 + CCA**
- 蛋白-代谢物相关性矩阵，典型相关分析
- **通路联合富集**
- 蛋白KEGG富集，代谢物通路富集，共同通路识别
- **因果调控网络**
- 蛋白→代谢物调控，反馈回路识别
- **多组学Biomarker发现**
- 联合评分，ROC分析准备

### Framework 5: 三组学修饰串扰
- **修饰位点-蛋白表达关联分类**
- Triple-Active/Suppressed，Phospho/Lact-Specific，Modification-Coordinated/Opposite
- **磷酸化-乳酸化位点级串扰**
- 位置接近度分析，共调控位点对
- **酶-底物网络构建**
- 激酶-底物，乳酸化酶-底物
- **激酶/乳酸化酶活性推断**
- 底物富集评分，活性状态判定

## 注意事项

- **数据格式**: 所有输入文件第一列应为基因/蛋白ID，列名应为样本名
- **缺失值**: 分析脚本会自动跳过缺失值，但建议先进行数据质控
- **物种**: 默认设置为human，如需其他物种请修改 `config/pipeline_config.yaml`
- **内存**: 大规模数据建议64GB以上内存

## 参考

- Subramanian et al. (2005) Gene set enrichment analysis. PNAS
- KSEA: Multiple kinase substrate analysis. Sci Signal
- Lactylation: Dai et al. (2021) Lactylation of PKM2... Nature
- STTT, NC, JECCR - 参考客户发表文献

## 联系方式

- 课题设计咨询：陈帅通团队
- 生信分析流程：OpenClaw AI助手
