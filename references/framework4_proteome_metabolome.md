# 框架4：蛋白质组 + 代谢组 联合分析

## 背景

代谢组是表型层面的直接反映，蛋白是功能执行者。两者联合是连接基因型与表型的桥梁。

**核心价值：**
- 揭示代谢重编程的分子机制
- 发现疾病biomarker（多组学panel）
- 药靶发现（代谢酶）
- 临床转化（代谢标志物）

## 分析流程

### Step 1: 数据独立分析

```python
# 蛋白组
proteome_df = pd.read_csv('proteome.tsv', sep='\t')
dep_df = proteome_df[proteome_df['P-value'] < 0.05]  # 差异蛋白

# 代谢组
metabolome_df = pd.read_csv('metabolome.tsv', sep='\t')
dem_df = metabolome_df[metabolome_df['P-value'] < 0.05]  # 差异代谢物
```

### Step 2: Spearman相关性分析

```python
import seaborn as sns
from scipy.stats import spearmanr

def correlation_analysis(proteome_df, metabolome_df, min_corr=0.7, pval_thresh=0.05):
    """
    蛋白-代谢物Spearman相关性分析
    """
    # 获取差异蛋白和差异代谢物
    dep_genes = dep_df['Gene'].unique()
    dem_mets = dem_df['Metabolite'].unique()
    
    # 构建表达矩阵
    prot_expr = proteome_df[proteome_df['Gene'].isin(dep_genes)].set_index('Gene')
    met_expr = metabolome_df[metabolome_df['Metabolite'].isin(dem_mets)].set_index('Metabolite')
    
    # 计算相关性矩阵
    corr_matrix = pd.DataFrame(index=prot_expr.index, columns=met_expr.index)
    pval_matrix = pd.DataFrame(index=prot_expr.index, columns=met_expr.index)
    
    for prot in prot_expr.index:
        for met in met_expr.index:
            rho, pval = spearmanr(prot_expr.loc[prot], met_expr.loc[met])
            corr_matrix.loc[prot, met] = rho
            pval_matrix.loc[prot, met] = pval
    
    # 筛选显著相关对
    significant_pairs = []
    for prot in corr_matrix.index:
        for met in corr_matrix.columns:
            rho = float(corr_matrix.loc[prot, met])
            pval = float(pval_matrix.loc[prot, met])
            if abs(rho) > min_corr and pval < pval_thresh:
                significant_pairs.append({
                    'Protein': prot,
                    'Metabolite': met,
                    'Correlation': rho,
                    'P-value': pval
                })
    
    return pd.DataFrame(significant_pairs).sort_values('Correlation', key=abs, ascending=False)
```

### Step 3: 典型相关分析（CCA）

```python
from sklearn.cross_decomposition import CCA
import numpy as np

def cca_analysis(proteome_matrix, metabolome_matrix, n_components=2):
    """
    典型相关分析：找出蛋白组和代谢组之间的最大相关性组合
    """
    # 标准化
    prot_std = (proteome_matrix - proteome_matrix.mean()) / proteome_matrix.std()
    met_std = (metabolome_matrix - metabolome_matrix.mean()) / metabolome_matrix.std()
    
    # CCA
    cca = CCA(n_components=n_components)
    prot_cc, met_cc = cca.fit_transform(prot_std, met_std)
    
    # 结果
    print(f"CCA correlations: {cca.coef_.shape}")
    print(f"Canonical correlations: {cca.score(prot_std, met_std):.3f}")
    
    return prot_cc, met_cc, cca

def visualize_cca(prot_cc, met_cc, sample_groups):
    """
    可视化CCA结果
    """
    plt.figure(figsize=(8, 6))
    
    # 样本在CCA空间的分布
    for group in set(sample_groups):
        mask = np.array(sample_groups) == group
        plt.scatter(prot_cc[mask, 0], prot_cc[mask, 1], label=f'{group} (Protein)', alpha=0.5)
        plt.scatter(met_cc[mask, 0], met_cc[mask, 1], marker='x', s=100, label=f'{group} (Metabolite)')
    
    plt.xlabel('CCA1')
    plt.ylabel('CCA2')
    plt.legend()
    plt.title('Canonical Correlation Analysis')
    plt.savefig('cca_plot.png', dpi=300)
```

### Step 4: 通路联合富集

```python
def pathway_integration_analysis(dep_genes, dem_mets, kegg_pathways):
    """
    同时对蛋白和代谢物进行通路富集，寻找共同通路
    """
    from scipy.stats import hypergeometric
    
    results = []
    for pathway_name, pathway_genes in kegg_pathways.items():
        # 获取该通路中的代谢物（如果有）
        pathway_mets = dem_mets[dem_mets.isin(pathway_genes)] if 'dem_mets' in dir() else []
        
        # 蛋白富集
        overlap_proteins = set(dep_genes) & set(pathway_genes)
        
        if len(overlap_proteins) > 0:
            results.append({
                'Pathway': pathway_name,
                'Type': 'Protein+Metabolite' if len(pathway_mets) > 0 else 'Protein only',
                'Protein_overlap': len(overlap_proteins),
                'Metabolite_overlap': len(pathway_mets),
                'Genes': list(overlap_proteins)
            })
    
    return pd.DataFrame(results)
```

### Step 5: 因果调控网络

```python
import networkx as nx

def build_regulatory_network(corr_pairs, enzyme_metabolite_db):
    """
    构建蛋白-代谢物调控网络
    - 酶-底物关系（已知）
    - 相关性关系（候选调控）
    """
    G = nx.DiGraph()
    
    # 添加酶-底物边（已知调控）
    for enzyme, metabolites in enzyme_metabolite_db.items():
        for met in metabolites:
            G.add_edge(enzyme, met, relation='enzyme_substrate', color='blue')
    
    # 添加显著相关边（候选调控）
    for _, row in corr_pairs.iterrows():
        # 检查是否已经是酶-底物关系
        if not G.has_edge(row['Protein'], row['Metabolite']):
            G.add_edge(
                row['Protein'], 
                row['Metabolite'], 
                relation='correlation',
                weight=abs(row['Correlation']),
                color='gray'
            )
    
    return G
```

### Step 6: 多组学Biomarker发现

```python
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LassoCV
from sklearn.model_selection import cross_val_score
import numpy as np

def biomarker_discovery(proteome_df, metabolome_df, labels, n_features=20):
    """
    多组学联合biomarker发现
    使用随机森林或LASSO回归
    """
    # 合并数据
    X = proteome_df.merge(metabolome_df, left_index=True, right_index=True)
    y = labels
    
    # 特征选择：LASSO
    lasso = LassoCV(cv=5)
    lasso.fit(X, y)
    
    # 获取重要特征
    feature_importance = pd.DataFrame({
        'Feature': X.columns,
        'Coefficient': lasso.coef_
    }).sort_values('Coefficient', key=abs, ascending=False)
    
    # 随机森林验证
    top_features = feature_importance.head(n_features)['Feature'].tolist()
    rf = RandomForestClassifier(n_estimators=100)
    scores = cross_val_score(rf, X[top_features], y, cv=5)
    
    print(f"Random Forest CV Accuracy: {scores.mean():.3f} ± {scores.std():.3f}")
    
    return feature_importance, top_features
```

## 代谢组数据库

| 数据库 | 描述 |
|--------|------|
| KEGG | 代谢通路 |
| HMDB | 人体代谢组 |
| MetaboAnalyst | 代谢分析工具 |
| MetID | 代谢物鉴定 |

## 典型应用

| 场景 | 研究问题 | 方法 |
|------|---------|------|
| 肿瘤代谢 | 糖酵解 vs 氧化磷酸化切换 | PCC/CCA |
| 糖尿病 | 胰岛素抵抗机制 | 相关性网络 |
| 药理 | 药物作用代谢标志物 | Biomarker panel |
| 营养 | 饮食干预代谢响应 | 时间序列分析 |

## 可视化

```python
def visualize_metabolome_proteome(corr_pairs):
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))
    
    # 1. 相关性热图
    ax1 = axes[0, 0]
    top_proteins = corr_pairs.head(30)['Protein'].unique()
    top_mets = corr_pairs.head(30)['Metabolite'].unique()
    corr_subset = corr_matrix.loc[top_proteins, top_mets].astype(float)
    sns.heatmap(corr_subset, cmap='RdBu_r', center=0, ax=ax1, annot=False)
    ax1.set_title('Protein-Metabolite Correlation Heatmap')
    
    # 2. 相关性分布
    ax2 = axes[0, 1]
    ax2.hist(corr_pairs['Correlation'], bins=30, edgecolor='black')
    ax2.set_xlabel('Correlation')
    ax2.set_ylabel('Count')
    ax2.set_title('Distribution of Correlations')
    
    # 3. Top相关性
    ax3 = axes[1, 0]
    top10 = corr_pairs.head(10)
    colors = ['red' if x > 0 else 'blue' for x in top10['Correlation']]
    ax3.barh(range(len(top10)), top10['Correlation'], color=colors)
    ax3.set_yticks(range(len(top10)))
    ax3.set_yticklabels([f"{p}-{m}" for p, m in zip(top10['Protein'], top10['Metabolite'])])
    ax3.set_xlabel('Correlation')
    ax3.set_title('Top 10 Protein-Metabolite Correlations')
    
    # 4. 网络图
    ax4 = axes[1, 1]
    G = build_regulatory_network(corr_pairs.head(50), {})
    pos = nx.spring_layout(G, k=2)
    nx.draw_networkx(G, pos, ax=ax4, node_size=100, font_size=5, alpha=0.7)
    ax4.set_title('Regulatory Network')
    
    plt.tight_layout()
    plt.savefig('metabolome_proteome_analysis.png', dpi=300)
```

## 参考文献

- Zamboni et al., Anal Chem 2015 (Metadata analysis)
- Hasin et al., Genome Biol 2017 (Multi-omics metabolomics)
- Mardinoglu et al., Cell 2017 (Integration of proteomics and metabolomics)
