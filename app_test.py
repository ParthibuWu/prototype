# app.py
import streamlit as st
import pandas as pd

from plots import (
    plot_gc_distribution,
    plot_length_distribution,
    plot_base_composition,
    plot_quality_distribution,
)

from processors import (
    get_sequence_stats,
    filter_fastq,
    process_sequences_universal,
)

st.set_page_config(page_title="Bioinformatics Sequence QC", layout="wide")
st.title("GenesisQC")

# ---- Sidebar: controls ----
st.sidebar.header("⚙️ Processing Options")

processing_mode = st.sidebar.selectbox(
    "Processing Mode",
    ["Statistics Only", "FASTQ Filtering", "Full Processing"]
)

denom_mode = st.sidebar.selectbox(
    "GC Denominator Mode",
    ["canonical", "raw"],
    help="canonical: denom = count(ATGC); raw: denom = total length including N/ambiguity"
)

max_records = st.sidebar.number_input(
    "Max records to load (Full/Filter mode)",
    min_value=10, max_value=5000, value=200, step=50
)

uploaded_file = st.file_uploader(
    "Upload FASTA/FASTQ/GenBank/EMBL (optionally .gz/.bz2)",
    type=["fa", "fasta", "fq", "fastq", "gb", "gbk", "embl", "gz", "bz2"]
)

if uploaded_file is None:
    st.info("Upload a file to begin.")
    st.stop()

# ---- File info metrics ----
col1, col2, col3 = st.columns(3)
with col1:
    st.metric("Filename", uploaded_file.name)
with col2:
    st.metric("File Size", f"{uploaded_file.size / 1024:.2f} KB")
with col3:
    st.metric("Type", uploaded_file.type or "Unknown")

st.markdown("---")

content = uploaded_file.getvalue()

# ---- Shared helper: render Phred summary metrics ----
def _render_phred_metrics(sequences: list) -> None:
    """Render a Phred-score summary row if quality fields are present."""
    if not sequences or "Avg_quality" not in sequences[0]:
        return

    avg_q_vals   = [s["Avg_quality"]   for s in sequences if s.get("Avg_quality")   is not None]
    min_q_vals   = [s["Min_quality"]   for s in sequences if s.get("Min_quality")   is not None]
    q20_pct_vals = [s["Q20_percent"]   for s in sequences if s.get("Q20_percent")   is not None]
    q30_pct_vals = [s["Q30_percent"]   for s in sequences if s.get("Q30_percent")   is not None]

    if not avg_q_vals:
        return

    st.subheader("📈 Phred Quality Summary")
    qc1, qc2, qc3, qc4 = st.columns(4)
    with qc1:
        st.metric("Mean Avg Quality", f"{sum(avg_q_vals)/len(avg_q_vals):.2f}")
    with qc2:
        st.metric("Lowest Min Quality", f"{min(min_q_vals):.0f}")
    with qc3:
        st.metric("Mean Q20%", f"{sum(q20_pct_vals)/len(q20_pct_vals):.1f}%")
    with qc4:
        st.metric("Mean Q30%", f"{sum(q30_pct_vals)/len(q30_pct_vals):.1f}%")


# =========================
# Mode 1: Statistics Only
# =========================
if processing_mode == "Statistics Only":
    if st.button("📊 Calculate Statistics", type="primary"):
        with st.spinner("Processing..."):
            result = get_sequence_stats(content, uploaded_file.name, mode=denom_mode)

        st.success("✅ Statistics calculated successfully!")

        is_fastq = result["format"] == "fastq"

        # Core metrics
        n_cols = 5 if is_fastq else 4
        cols = st.columns(n_cols)
        with cols[0]:
            st.metric("Total Sequences", result["total_sequences"])
        with cols[1]:
            st.metric("Total Bases", f"{result['total_bases']:,}")
        with cols[2]:
            st.metric("Avg Length", f"{result['average_length']:.2f} bp")
        with cols[3]:
            st.metric("Avg GC%", f"{result['average_gc_content']:.2f}%")
        if is_fastq and "average_quality" in result:
            with cols[4]:
                st.metric("Avg Phred Quality", f"{result['average_quality']:.2f}")

        st.info(f"📁 Format: {result['format'].upper()} | Compression: {result['compression'].upper()}")

        if is_fastq and "average_quality" in result:
            st.markdown("---")
            st.subheader("📋 Quality Score Interpretation")
            avg_q = result["average_quality"]
            if avg_q >= 30:
                st.success(f"Q{avg_q:.0f} — Excellent quality (error rate < 0.1%)")
            elif avg_q >= 20:
                st.warning(f"Q{avg_q:.0f} — Acceptable quality (error rate < 1%)")
            else:
                st.error(f"Q{avg_q:.0f} — Poor quality (error rate > 1%) — consider trimming")


# =========================
# Mode 2: FASTQ Filtering
# =========================
elif processing_mode == "FASTQ Filtering":
    filter_ids = st.text_area(
        "Filter by Sequence IDs (comma-separated)",
        placeholder="seq1, seq2, seq3",
        help="Leave empty to process all sequences."
    )

    if st.button("🔍 Process & Filter", type="primary"):
        with st.spinner("Processing..."):
            result = filter_fastq(
                content, uploaded_file.name,
                filter_ids or None,
                mode=denom_mode,
                max_records=int(max_records)
            )

        st.success("✅ Processing complete!")

        c1, c2, c3 = st.columns(3)
        with c1:
            st.metric("Total Sequences", result["total_sequences"])
        with c2:
            st.metric("Filtered", "Yes" if result["filtered"] else "No")
        with c3:
            st.metric("Filter Count", result.get("filter_count", 0))

        sequences = result.get("sequences", [])

        if sequences:
            # Phred summary
            _render_phred_metrics(sequences)

            st.markdown("---")
            df = pd.DataFrame(sequences)
            st.subheader("📄 Filtered Sequences")
            st.dataframe(df, use_container_width=True, height=400)

            csv = df.to_csv(index=False).encode("utf-8")
            st.download_button(
                label="⬇️ Download Results (CSV)",
                data=csv,
                file_name="filtered_sequences.csv",
                mime="text/csv",
            )

            st.subheader("📊 Visualizations")
            v1, v2 = st.columns(2)
            with v1:
                fig = plot_gc_distribution(sequences)
                if fig:
                    st.plotly_chart(fig, use_container_width=True)
            with v2:
                fig = plot_length_distribution(sequences)
                if fig:
                    st.plotly_chart(fig, use_container_width=True)

            qfig = plot_quality_distribution(sequences)
            if qfig is not None:
                st.subheader("Quality Scores Distribution")
                st.plotly_chart(qfig, use_container_width=True)


# =========================
# Mode 3: Full Processing
# =========================
else:
    if st.button("🚀 Process File", type="primary"):
        with st.spinner("Processing sequences..."):
            result = process_sequences_universal(
                content, uploaded_file.name,
                mode=denom_mode,
                max_records=int(max_records)
            )

        st.success("✅ Done!")

        c1, c2, c3, c4 = st.columns(4)
        with c1:
            st.metric("Format", result["format"].upper())
        with c2:
            st.metric("Compression", result["compression"].upper())
        with c3:
            st.metric("Sequences", result.get("total_sequences", len(result.get("sequences", []))))
        with c4:
            st.metric("Total Bases", f"{result.get('total_bases', 0):,}")

        sequences = result.get("sequences", [])
        if not sequences:
            st.warning("No sequences parsed.")
            st.stop()

        df = pd.DataFrame(sequences)

        tab1, tab2, tab3, tab4 = st.tabs(["📄 Data Table", "📊 Visualizations", "🔎 Details", "⬇️ Export"])

        # ---- Tab 1: Data Table ----
        with tab1:
            st.subheader("Summary Statistics")

            gc_col = next((c for c in ["GC_percent", "GC_content"] if c in df.columns), None)
            is_fastq_data = "Avg_quality" in df.columns

            # Base metrics
            base_cols = st.columns(5 if is_fastq_data else 3)
            with base_cols[0]:
                st.metric("Mean Length", f"{df['Length'].mean():.2f} bp" if "Length" in df else "N/A")
            with base_cols[1]:
                st.metric("Mean GC%", f"{df[gc_col].mean():.2f}%" if gc_col else "N/A")
            with base_cols[2]:
                st.metric("Std Dev GC%", f"{df[gc_col].std():.2f}%" if gc_col else "N/A")

            if is_fastq_data:
                with base_cols[3]:
                    st.metric("Mean Avg Quality", f"{df['Avg_quality'].mean():.2f}")
                with base_cols[4]:
                    st.metric("Mean Q30%", f"{df['Q30_percent'].mean():.1f}%")

            st.dataframe(df, use_container_width=True, height=500)

            # Phred per-record summary table (FASTQ only)
            if is_fastq_data:
                st.markdown("---")
                _render_phred_metrics(sequences)

                phred_cols = ["ID", "Length", "Avg_quality", "Min_quality", "Max_quality",
                              "Q20_bases", "Q30_bases", "Q20_percent", "Q30_percent"]
                available = [c for c in phred_cols if c in df.columns]
                st.subheader("🔬 Per-Record Phred Stats")
                st.dataframe(
                    df[available].style.format({
                        "Avg_quality": "{:.2f}",
                        "Q20_percent": "{:.1f}",
                        "Q30_percent": "{:.1f}",
                    }),
                    use_container_width=True,
                    height=350,
                )

        # ---- Tab 2: Visualizations ----
        with tab2:
            left, right = st.columns(2)
            with left:
                fig = plot_gc_distribution(sequences)
                if fig:
                    st.plotly_chart(fig, use_container_width=True)
            with right:
                fig = plot_length_distribution(sequences)
                if fig:
                    st.plotly_chart(fig, use_container_width=True)

            if sequences and "A" in sequences[0] and "denom_used" in sequences[0]:
                fig = plot_base_composition(sequences)
                if fig:
                    st.plotly_chart(fig, use_container_width=True)

            qfig = plot_quality_distribution(sequences)
            if qfig is not None:
                st.subheader("📉 Quality Score Distribution")
                st.plotly_chart(qfig, use_container_width=True)

        # ---- Tab 3: Details ----
        with tab3:
            st.subheader("Detailed Sequence Information")

            id_field = "ID" if "ID" in df.columns else df.columns[0]
            idx = st.selectbox(
                "Select sequence to view",
                options=list(range(len(sequences))),
                format_func=lambda i: sequences[i].get(id_field, f"seq_{i}")
            )

            selected = sequences[idx]
            colA, colB = st.columns(2)

            with colA:
                st.markdown(f"**ID:** {selected.get('ID', 'N/A')}")
                st.markdown(f"**Description:** {selected.get('Description', selected.get('Title', 'N/A'))}")
                st.markdown(f"**Length:** {selected.get('Length', 'N/A')} bp")
                gc_val = selected.get("GC_percent", selected.get("GC_content", None))
                st.markdown(f"**GC%:** {f'{gc_val:.2f}' if gc_val is not None else 'N/A'}")

                # Phred stats block
                if "Avg_quality" in selected:
                    st.markdown("---")
                    st.markdown("**Phred Quality Stats**")
                    avg_q = selected["Avg_quality"]
                    st.markdown(f"- Avg Quality: **{avg_q:.2f}**")
                    st.markdown(f"- Min Quality: **{selected.get('Min_quality', 'N/A')}**")
                    st.markdown(f"- Max Quality: **{selected.get('Max_quality', 'N/A')}**")
                    st.markdown(f"- Q20 bases: **{selected.get('Q20_bases', 'N/A')}** ({selected.get('Q20_percent', 0):.1f}%)")
                    st.markdown(f"- Q30 bases: **{selected.get('Q30_bases', 'N/A')}** ({selected.get('Q30_percent', 0):.1f}%)")

                    # Quality badge
                    if avg_q >= 30:
                        st.success(f"Q{avg_q:.0f} — Excellent (error rate < 0.1%)")
                    elif avg_q >= 20:
                        st.warning(f"Q{avg_q:.0f} — Acceptable (error rate < 1%)")
                    else:
                        st.error(f"Q{avg_q:.0f} — Poor (error rate > 1%)")

            with colB:
                seq_str = selected.get("Sequence", "")
                if seq_str:
                    if len(seq_str) > 500:
                        st.text_area("Sequence (first 500 bp)", seq_str[:500] + "...", height=200)
                        st.caption(f"Full sequence length: {len(seq_str)} bp")
                    else:
                        st.text_area("Sequence", seq_str, height=200)
                else:
                    st.info("Sequence string not stored in this mode (memory-saving).")

        # ---- Tab 4: Export ----
        with tab4:
            st.subheader("Export Options")

            csv = df.to_csv(index=False).encode("utf-8")
            st.download_button(
                label="⬇️ Download Full Table (CSV)",
                data=csv,
                file_name=f"{uploaded_file.name}_processed.csv",
                mime="text/csv",
            )

            # Separate Phred-only export for FASTQ
            if "Avg_quality" in df.columns:
                phred_cols = ["ID", "Length", "Avg_quality", "Min_quality", "Max_quality",
                              "Q20_bases", "Q30_bases", "Q20_percent", "Q30_percent"]
                available = [c for c in phred_cols if c in df.columns]
                phred_csv = df[available].to_csv(index=False).encode("utf-8")
                st.download_button(
                    label="⬇️ Download Phred Stats Only (CSV)",
                    data=phred_csv,
                    file_name=f"{uploaded_file.name}_phred_stats.csv",
                    mime="text/csv",
                )