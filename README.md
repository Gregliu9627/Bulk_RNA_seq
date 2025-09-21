# Bulk RNA-seq Analysis Portal

该项目提供一个运行在服务器上的 Streamlit 网页，用于完成从原始 FASTQ 到下游差异分析、可视化和通路注释的整个 bulk RNA-seq 流程。

## 功能概览

- **原始数据处理**：基于 Salmon 的定量流程，可选 fastp 质控，自动汇总为 gene-level counts。
- **差异表达分析**：通过内置的 `scripts/run_deseq2.R` 调用 DESeq2，生成标准差异基因列表及规范化矩阵。
- **PCA 可视化**：支持按方差选择基因、标准化、元数据着色。
- **火山图**：可调节 log2FC/padj 阈值、颜色、x 轴范围等。
- **热图与聚类**：基于方差筛选基因，支持自定义色板。
- **富集分析**：调用 gseapy 进行 prerank GSEA，基因集文件放置于 `data/reference/gene_sets/`。
- **趋势分析**：针对元数据中的数值列进行线性趋势检测。

## 文件结构

```
.
├── app/                # Streamlit 前端
├── bulk_rna_seq/       # Python 工作流与分析模块
├── data/
│   ├── raw/            # 放置原始 FASTQ 及 samplesheet.csv、metadata.csv
│   ├── processed/      # 流程输出的 counts 目录
│   └── reference/      # 参考资源（Salmon index, GTF, gene sets）
├── results/            # 下游分析输出
├── scripts/run_deseq2.R# 调用 DESeq2 的 R 脚本
└── requirements.txt
```

## 准备数据

1. 在 `data/raw/<dataset>` 中放置 FASTQ 文件，并创建 `samplesheet.csv`，格式示例：

```csv
sample,fastq_1,fastq_2
sample1,sample1_R1.fastq.gz,sample1_R2.fastq.gz
sample2,sample2_R1.fastq.gz,sample2_R2.fastq.gz
```

2. 在同一目录下提供 `metadata.csv`，第一列为样本 ID，其余列为实验信息（如条件、批次、时间等）。
3. 将 Salmon index 放在 `data/reference/salmon_index`，GTF 注释（可选）放在 `data/reference/annotation.gtf`。
4. 将基因集 `.gmt/.gmx` 文件放入 `data/reference/gene_sets/` 以用于 GSEA。

## 运行 Web 界面

```bash
pip install -r requirements.txt
streamlit run app/app.py --server.port 8501 --server.address 0.0.0.0
```

访问 `http://<服务器IP>:8501` 即可。

## 常见问题

- **缺少外部工具**：原始数据处理依赖 `fastp` 和 `salmon`，请确保二者在 `PATH` 中。
- **DESeq2 报错**：确认服务器上已安装 R 以及 DESeq2、optparse 等依赖。
- **gseapy 未安装**：若不需要 GSEA，可从 `requirements.txt` 中移除或忽略相关模块。
