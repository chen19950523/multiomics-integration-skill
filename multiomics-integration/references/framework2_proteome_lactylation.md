# 框架2：蛋白质组 + 乳酸化（Lactylation）联合分析

## 背景

乳酸化（Lactylation, Kla）是2019年发现的新型蛋白质翻译后修饰（PTM），由乳酸直接修饰赖氨酸残基，参与代谢重编程和表观遗传调控。

**核心研究价值：**
- 代谢-表观遗传桥梁
- 线粒体功能调控
- 肿瘤免疫代谢
- 胚胎发育

## 分析流程

### Step 1: 乳酸化蛋白组分析

**数据来源：** 抗Kla抗体富集 + LC-MS/MS
**常用软件：** MaxQuant (with Gla标签定义), pFind

```python
# 乳酸化位点数据格式
lactyl_df = pd.DataFrame({
    'ProteinID': ['P12345', 'P12345'],
    'Gene': ['GAPDH', 'GAPDH'],
    'Site': ['K184', 'K334'],
    'Position': [184, 334],
    'log2FC': [1.2, -0.8],
    'P-value': [0.001, 0.05],
    'Regulation': ['Up', 'Down']
})
```

### Step 2: 全局蛋白组数据关联

```python
def lactylation_protein_association(lactyl_df, proteome_df):
    """
    乳酸化位点与全局蛋白表达的关联分析
    """
    merged = lactyl_df.merge(
        proteome_df[['ProteinID', 'log2FC_protein', 'P-value_protein', 'Regulation_protein']], 
        on='ProteinID', 
        how='left'
    )
    return merged
```

### Step 3: 7种表达-修饰关系分类

| Type | 蛋白表达 | 乳酸化 | 生物学含义 |
|------|---------|--------|-----------|
| A | ↑ | ↑ | 表达升高伴随乳酸化增加（代谢感知）|
| B | ↑ | ↓ | 乳酸化被稀释（蛋白↑但修饰比例↓）|
| C | ↓ | ↑ | 乳酸化保护蛋白不被降解 |
| D | ↓ | ↓ | 蛋白减少，乳酸化也减少 |
| E | 不变 | ↑ | **纯乳酸化调控（功能激活）**|
| F | 不变 | ↓ | 乳酸化抑制 |
| G | 不变 | 不变 | 不受影响 |

```python
def classify_lactylation_types(merged_df, fc_thresh=0.58, p_thresh=0.05):
    """
    对蛋白-乳酸化关系进行7分类
    """
    conditions = {
        'A': (merged_df['Protein_FC'] > fc_thresh) & (merged_df['Lactyl_FC'] > fc_thresh),
        'B': (merged_df['Protein_FC'] > fc_thresh) & (merged_df['Lactyl_FC'] < -fc_thresh),
        'C': (merged_df['Protein_FC'] < -fc_thresh) & (merged_df['Lactyl_FC'] > fc_thresh),
        'D': (merged_df['Protein_FC'] < -fc_thresh) & (merged_df['Lactyl_FC'] < -fc_thresh),
        'E': (abs(merged_df['Protein_FC']) <= fc_thresh) & (merged_df['Lactyl_FC'] > fc_thresh),
        'F': (abs(merged_df['Protein_FC']) <= fc_thresh) & (merged_df['Lactyl_FC'] < -fc_thresh),
        'G': (abs(merged_df['Protein_FC']) <= fc_thresh) & (abs(merged_df['Lactyl_FC']) <= fc_thresh),
    }
    
    for type_name, mask in conditions.items():
        genes = merged_df[mask]['Gene'].tolist()
        print(f"Type {type_name}: {len(genes)} proteins")
    
    return conditions
```

### Step 4: 乳酸化Motif分析

```python
def lactylation_motif_analysis(site_sequences):
    """
    分析乳酸化位点周围的序列motif
    常用工具：MoMo (MEME Suite)
    """
    from collections import Counter
    
    # 提取-6到+6位置的氨基酸
    motifs = []
    for seq in site_sequences:
        if len(seq) >= 13:  # K在第7位
            motif = seq[6:13]  # K及其右侧
            motifs.append(motif)
    
    # 统计频率
    freq = Counter(motifs)
    return freq

# 常见乳酸化motif示例
# [K][K/P][K] - 偏好赖氨酸和脯氨酸
```

### Step 5: 乳酸化酶活性预测

```python
# Writer酶（写入）：HATs (KAT2A/KAT2B), CBP
# Eraser酶（擦除）：HDAC1-3, SIRT1-3

LACTYLATION_ENZYMES = {
    'Writer': ['KAT2A', 'KAT2B', 'CREBBP', 'EP300'],
    'Eraser': ['HDAC1', 'HDAC2', 'HDAC3', 'SIRT1', 'SIRT2', 'SIRT3']
}

def predict_lactylase_activity(target_proteins, enzyme_list):
    """
    预测哪些底物蛋白受特定乳酸化酶调控
    基于蛋白-蛋白相互作用预测
    """
    predicted_targets = {}
    for enzyme in enzyme_list:
        # 简化：检查酶与底物是否共表达
        predicted_targets[enzyme] = [
            p for p in target_proteins 
            if enzyme.lower() in p.lower() or p.lower() in enzyme.lower()
        ]
    return predicted_targets
```

### Step 6: 功能富集（重点通路）

```python
def lactylation_pathway_enrichment(lactylated_genes):
    """
    乳酸化蛋白的功能富集
    重点关注代谢和线粒体相关通路
    """
    pathways = {
        'Mitochondrial_function': ['MT-ATP6', 'MT-CO1', 'MT-CO2', 'NDUFA9', 'NDUFB8'],
        'Glycolysis': ['HK1', 'HK2', 'PGK1', 'PKM', 'LDHA'],
        'Gluconeogenesis': ['PCK1', 'PCK2', 'FBP1', 'G6PC'],
        'Fatty_acid_oxidation': ['ACADM', 'ACADL', 'CPT1A', 'CPT2'],
        'TCA_cycle': ['IDH1', 'IDH2', 'SDHA', 'SUCLA2'],
        'Histone_regulation': ['H1', 'H2A', 'H2B', 'H3', 'H4'],
    }
    
    results = []
    for pathway, genes in pathways.items():
        overlap = set(lactylated_genes) & set(genes)
        if overlap:
            results.append({
                'Pathway': pathway,
                'Overlap_count': len(overlap),
                'Genes': list(overlap)
            })
    
    return pd.DataFrame(results)
```

### Step 7: 位点功能预测

```python
def predict_lactylation_function(protein, site, domains):
    """
    预测乳酸化位点的功能影响
    - 是否在活性位点？
    - 是否在蛋白互作界面？
    - 是否在结构域内？
    """
    predictions = []
    for domain in domains:
        if domain['start'] <= site <= domain['end']:
            predictions.append(f"Site {site} is in {domain['name']} domain")
    
    # 在线工具：MusiteDeep, LabCaMP
    return predictions
```

## 可视化

```python
import seaborn as sns
import matplotlib.pyplot as plt

def visualize_lactylation_results(merged_df):
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))
    
    # 1. 蛋白表达 vs 乳酸化散点图
    ax1 = axes[0, 0]
    sns.scatterplot(data=merged_df, x='Protein_log2FC', y='Lactyl_log2FC', ax=ax1, alpha=0.5)
    ax1.axhline(0, color='red', linestyle='--', alpha=0.3)
    ax1.axvline(0, color='red', linestyle='--', alpha=0.3)
    ax1.set_title('Protein Expression vs Lactylation')
    
    # 2. Type分布饼图
    ax2 = axes[0, 1]
    type_counts = merged_df['Type'].value_counts()
    ax2.pie(type_counts.values, labels=type_counts.index, autopct='%1.1f%%')
    ax2.set_title('Distribution of Expression-Lactylation Types')
    
    # 3. 通路富集条形图
    ax3 = axes[1, 0]
    pathway_data = lactylation_pathway_enrichment(merged_df['Gene'].tolist())
    if not pathway_data.empty:
        sns.barplot(data=pathway_data, x='Overlap_count', y='Pathway', ax=ax3)
    ax3.set_title('Pathway Enrichment')
    
    # 4. 位点分布图（蛋白结构示意）
    ax4 = axes[1, 1]
    ax4.scatter(merged_df['Position'], [1]*len(merged_df), c=merged_df['Lactyl_log2FC'], cmap='RdBu_r')
    ax4.set_xlabel('Amino Acid Position')
    ax4.set_title('Lactylation Site Distribution')
    
    plt.tight_layout()
    plt.savefig('lactylation_analysis.png', dpi=300)
```

## 典型应用场景

| 疾病/模型 | 研究重点 | 预期发现 |
|----------|---------|---------|
| 肿瘤代谢 | Warburg effect + 乳酸化 | 糖酵解通路乳酸化异常 |
| 免疫细胞 | T细胞活化和分化 | 免疫检查点乳酸化调控 |
| 神经退行 | 线粒体功能 | 脑组织乳酸化谱改变 |
| 心血管 | 心肌缺血再灌注 | 线粒体蛋白乳酸化保护作用 |

## 参考文献

- Zhang et al., Nature 2019 (Lactylation discovery)
- Wang et al., Cell Metab 2021 (Lactylation in tumor immunity)
- Yang et al., Nat Cell Biol 2022 (Histone lactylation and pluripotency)
