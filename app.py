"""Streamlit application orchestrating bulk RNA-seq analysis."""
from __future__ import annotations

from typing import Optional

import numpy as np
import pandas as pd
import streamlit as st

from bulk_rna_seq import config
from bulk_rna_seq.downstream import (
    DownstreamAnalysisError,
    compute_trends,
    generate_heatmap,
    generate_volcano_plot,
    load_counts,
    load_metadata,
    perform_pca,
    run_deseq2,
    run_gsea,
)
from bulk_rna_seq.raw_processing import (
    RawProcessingConfig,
    RawProcessingError,
    RawProcessingPipeline,
    list_available_datasets,
    list_processed_datasets,
)

st.set_page_config(page_title="Bulk RNA-seq Analysis Portal", layout="wide")


@st.cache_data(show_spinner=False)
def _load_counts(dataset: str) -> pd.DataFrame:
    return load_counts(dataset)


@st.cache_data(show_spinner=False)
def _load_metadata(dataset: str) -> pd.DataFrame:
    return load_metadata(dataset)


def _load_deseq_results(dataset: str) -> Optional[pd.DataFrame]:
    result_path = config.RESULTS_DIR / dataset / "deseq2" / "deseq2_results.csv"
    if result_path.exists():
        return pd.read_csv(result_path, index_col=0)
    return None


st.title("Bulk RNA-seq 分析门户")
st.caption("在服务器上一站式完成从FASTQ到生物学解释的所有步骤。")

raw_tab, downstream_tab = st.tabs(["原始数据处理", "下游分析"])


with raw_tab:
    st.header("原始数据处理：FASTQ → counts")
    datasets = list_available_datasets()
    if not datasets:
        st.info("请先将FASTQ数据和samplesheet.csv放入 data/raw/<dataset> 目录。")
    else:
        dataset = st.selectbox("选择数据集", datasets)
        threads = st.slider("线程数", min_value=1, max_value=32, value=config.DEFAULT_THREADS)
        run_fastp = st.checkbox("执行fastp质控和剪切", value=True)
        library_type = st.text_input("Salmon library type", value="A")
        force = st.checkbox("覆盖已有结果", value=False)

        if st.button("运行原始数据处理", type="primary"):
            params = RawProcessingConfig(
                dataset=dataset,
                threads=threads,
                run_fastp=run_fastp,
                library_type=library_type,
                force=force,
            )
            pipeline = RawProcessingPipeline(params)
            try:
                counts_path = pipeline.run()
            except RawProcessingError as exc:
                st.error(str(exc))
            else:
                st.success(f"完成！counts文件位于 {counts_path}")
                st.cache_data.clear()


with downstream_tab:
    st.header("下游分析模块")
    processed_datasets = list_processed_datasets()
    if not processed_datasets:
        st.warning("尚未检测到处理好的数据集，请先运行原始数据处理。")
    else:
        dataset = st.selectbox("选择已处理数据集", processed_datasets, key="processed_dataset")
        counts = _load_counts(dataset)
        metadata = _load_metadata(dataset)

        st.subheader("样本信息预览")
        st.dataframe(metadata)

        de_tab, pca_tab, heatmap_tab, gsea_tab, trend_tab = st.tabs(
            ["差异分析 (DESeq2)", "PCA", "聚类与热图", "富集分析", "趋势分析"]
        )

        with de_tab:
            st.markdown("### DESeq2 差异分析")
            design_formula = st.text_input(
                "Design formula", value=config.DEFAULT_DESIGN_FORMULA, key="design_formula"
            )
            condition_column = st.selectbox(
                "条件列 (用于contrast)", metadata.columns.tolist(), index=0
            )
            groups = metadata[condition_column].unique().tolist()
            control_group = st.selectbox("对照组", groups, index=0)
            treatment_group = st.selectbox("处理组", groups, index=min(1, len(groups) - 1))
            if st.button("运行DESeq2", key="run_deseq2"):
                try:
                    result = run_deseq2(
                        dataset=dataset,
                        counts=counts,
                        metadata=metadata,
                        design_formula=design_formula,
                        condition_column=condition_column,
                        control_group=control_group,
                        treatment_group=treatment_group,
                    )
                except DownstreamAnalysisError as exc:
                    st.error(str(exc))
                else:
                    st.success(f"DESeq2 完成，结果保存在 {result.result_file}")

            deseq_results = _load_deseq_results(dataset)
            if deseq_results is not None:
                st.markdown("#### 结果预览")
                st.dataframe(deseq_results.head(20))

                st.markdown("#### 火山图参数")
                col1, col2, col3 = st.columns(3)
                with col1:
                    log2fc_threshold = st.slider("log2FC 阈值", 0.0, 5.0, 1.0, 0.1)
                    padj_threshold = st.slider("padj 阈值", 0.0, 0.5, 0.05, 0.01)
                with col2:
                    max_log2fc = st.slider("log2FC 范围", 0.5, 10.0, 5.0, 0.5)
                    color_up = st.color_picker("上调颜色", "#d73027")
                with col3:
                    color_down = st.color_picker("下调颜色", "#4575b4")
                    color_ns = st.color_picker("非显著颜色", "#808080")

                fig = generate_volcano_plot(
                    deseq_results,
                    log2fc_threshold=log2fc_threshold,
                    padj_threshold=padj_threshold,
                    max_log2fc=max_log2fc,
                    color_up=color_up,
                    color_down=color_down,
                    color_ns=color_ns,
                )
                st.plotly_chart(fig, use_container_width=True)

                top_n = st.slider("导出Top基因数量", min_value=10, max_value=500, value=100, step=10)
                st.download_button(
                    "下载差异基因列表",
                    data=deseq_results.sort_values("padj").head(top_n).to_csv().encode("utf-8"),
                    file_name=f"{dataset}_deseq2_top{top_n}.csv",
                    mime="text/csv",
                )
            else:
                st.info("尚未运行DESeq2 或结果不存在。")

        with pca_tab:
            st.markdown("### PCA")
            top_n = st.slider("按方差选择前N个基因", min_value=50, max_value=2000, value=500, step=50)
            color_by = st.selectbox("颜色映射列", ["无"] + metadata.columns.tolist())
            hover_cols = st.multiselect("悬停显示字段", metadata.columns.tolist(), default=metadata.columns.tolist())
            scale = st.checkbox("标准化表达矩阵", value=True)
            if st.button("绘制PCA", key="run_pca"):
                fig = perform_pca(
                    counts=counts,
                    metadata=metadata,
                    top_n=top_n,
                    scale=scale,
                    color_by=None if color_by == "无" else color_by,
                    hover_data=hover_cols,
                )
                st.plotly_chart(fig, use_container_width=True)

        with heatmap_tab:
            st.markdown("### 基因聚类与热图")
            top_n = st.slider("热图包含的基因数量", min_value=20, max_value=200, value=50, step=10)
            color_scale = st.selectbox(
                "颜色方案",
                ["RdBu", "Viridis", "Plasma", "Inferno", "Magma", "Blues", "Greens"],
                index=0,
            )
            if st.button("生成热图", key="run_heatmap"):
                heatmap_fig = generate_heatmap(counts=counts, metadata=metadata, top_n=top_n, color_scale=color_scale)
                st.plotly_chart(heatmap_fig, use_container_width=True)

        with gsea_tab:
            st.markdown("### 富集分析")
            if _load_deseq_results(dataset) is None:
                st.info("请先运行DESeq2以生成排序列表。")
            else:
                deseq_results = _load_deseq_results(dataset)
                assert deseq_results is not None
                ranking_metric = st.selectbox(
                    "排名指标",
                    ["log2FoldChange", "stat", "log2FC * -log10(padj)"],
                )
                if ranking_metric == "log2FC * -log10(padj)":
                    padj = deseq_results["padj"].replace(0, pd.NA).astype(float)
                    ranking = deseq_results["log2FoldChange"] * -np.log10(padj)
                else:
                    ranking = deseq_results[ranking_metric]
                ranking = ranking.replace([np.inf, -np.inf], pd.NA).dropna()
                ranking = ranking.sort_values(ascending=False)

                gene_sets = [p.name for p in config.GENE_SET_DIR.glob("*") if p.is_file()]
                if not gene_sets:
                    st.warning("未找到基因集文件，请将.gmt或.gmx文件放入 data/reference/gene_sets/")
                else:
                    gene_set = st.selectbox("选择基因集文件", gene_sets)
                    min_size = st.slider("最小基因集大小", 5, 200, 15)
                    max_size = st.slider("最大基因集大小", 50, 2000, 500)
                    if st.button("运行GSEA", key="run_gsea"):
                        try:
                            outdir = run_gsea(
                                dataset=dataset,
                                ranking=ranking,
                                gene_set=gene_set,
                                min_size=min_size,
                                max_size=max_size,
                            )
                        except DownstreamAnalysisError as exc:
                            st.error(str(exc))
                        else:
                            st.success(f"GSEA 完成，结果保存在 {outdir}")

        with trend_tab:
            st.markdown("### 趋势分析")
            numeric_columns = metadata.select_dtypes(include=["number"]).columns.tolist()
            if not numeric_columns:
                st.info("元数据中没有数值型列，无法进行趋势分析。")
            else:
                numeric_column = st.selectbox("选择用于趋势分析的数值列", numeric_columns)
                top_n = st.slider("展示前N个基因", min_value=10, max_value=200, value=50, step=10)
                if st.button("计算趋势", key="run_trends"):
                    trends = compute_trends(counts=counts, metadata=metadata, numeric_column=numeric_column, top_n=top_n)
                    st.dataframe(trends)
                    st.download_button(
                        "下载趋势结果",
                        data=trends.to_csv(index=False).encode("utf-8"),
                        file_name=f"{dataset}_trends_{numeric_column}.csv",
                        mime="text/csv",
                    )
