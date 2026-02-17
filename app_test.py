# app.py
import streamlit as st
import pandas as pd

from plots import (
    plot_gc_distribution,
    plot_length_distribution,
    plot_base_composition,
    plot_quality_distribution,
)

# ---- IMPORT YOUR BACKEND FUNCTIONS ----
# These should already exist in your project (from your screenshots / earlier synthesis).
# If names differ, map them in the "Mapping" section at bottom.
from processors import (
    get_sequence_stats,
    filter_fastq,
    process_sequences_universal,
)

st.set_page_config(page_title="Bioinformatics Sequence QC", layout="wide")
st.title("🧬 Bioinformatics Sequence QC Dashboard")

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
    type=["fa","fasta","fq","fastq","gb","gbk","embl","gz","bz2"]
)

if uploaded_file is None:
    st.info("Upload a file to begin.")
    st.stop()

# ---- File info metrics (your screenshot feature) ----
col1, col2, col3 = st.columns(3)
with col1:
    st.metric("Filename", uploaded_file.name)
with col2:
    st.metric("File Size", f"{uploaded_file.size / 1024:.2f} KB")
with col3:
    st.metric("Type", uploaded_file.type or "Unknown")

st.markdown("---")

# We read bytes once; everything else should stream-parse internally
content = uploaded_file.getvalue()

# =========================
# Mode 1: Statistics Only
# =========================
if processing_mode == "Statistics Only":
    if st.button("📊 Calculate Statistics", type="primary"):
        with st.spinner("Processing..."):
            result = get_sequence_stats(content, uploaded_file.name, mode=denom_mode)

        st.success("✅ Statistics calculated successfully!")

        # Display stats metrics
        c1, c2, c3, c4 = st.columns(4)
        with c1:
            st.metric("Total Sequences", result["total_sequences"])
        with c2:
            st.metric("Total Bases", f"{result['total_bases']:,}")
        with c3:
            st.metric("Avg Length", f"{result['average_length']:.2f} bp")
        with c4:
            st.metric("Avg GC%", f"{result['average_gc_content']:.2f}%")

        st.info(f"📁 Format: {result['format'].upper()} | Compression: {result['compression'].upper()}")

# =========================
# Mode 2: FASTQ Filtering
# =========================
elif processing_mode == "FASTQ Filtering":
    filter_ids = st.text_area(
        "Filter by Sequence IDs (comma-separated)",
        placeholder="seq1, seq2, seq3",
        help="Leave empty to process all sequences (but this is meant for FASTQ)."
    )

    if st.button("🔍 Process & Filter", type="primary"):
        with st.spinner("Processing..."):
            result = filter_fastq(content, uploaded_file.name, filter_ids or None, mode=denom_mode, max_records=int(max_records))

        st.success("✅ Processing complete!")

        # Summary metrics
        c1, c2, c3 = st.columns(3)
        with c1:
            st.metric("Total Sequences", result["total_sequences"])
        with c2:
            st.metric("Filtered", "Yes" if result["filtered"] else "No")
        with c3:
            st.metric("Filter Count", result.get("filter_count", 0))

        if result["sequences"]:
            df = pd.DataFrame(result["sequences"])
            st.subheader("📄 Filtered Sequences")
            st.dataframe(df, use_container_width=True, height=400)

            csv = df.to_csv(index=False).encode("utf-8")
            st.download_button(
                label="⬇️ Download Results (CSV)",
                data=csv,
                file_name="filtered_sequences.csv",
                mime="text/csv",
            )

            # Visualizations (like your screenshot)
            st.subheader("📊 Visualizations")
            v1, v2 = st.columns(2)
            with v1:
                fig = plot_gc_distribution(result["sequences"])
                if fig: st.plotly_chart(fig, use_container_width=True)
            with v2:
                fig = plot_length_distribution(result["sequences"])
                if fig: st.plotly_chart(fig, use_container_width=True)

            # Quality plot if available
            qfig = plot_quality_distribution(result["sequences"])
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

        # Header metrics row (format/compression/sequences/bases) like your screenshot
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
            s1, s2, s3 = st.columns(3)
            with s1:
                st.metric("Mean Length", f"{df['Length'].mean():.2f} bp" if "Length" in df else "N/A")
            with s2:
                gc_col = "GC_percent" if "GC_percent" in df.columns else ("GC_content" if "GC_content" in df.columns else None)
                st.metric("Mean GC%", f"{df[gc_col].mean():.2f}%" if gc_col else "N/A")
            with s3:
                st.metric("Std Dev GC%", f"{df[gc_col].std():.2f}%" if gc_col else "N/A")

            st.dataframe(df, use_container_width=True, height=500)

        # ---- Tab 2: Visualizations ----
        with tab2:
            left, right = st.columns(2)
            with left:
                fig = plot_gc_distribution(sequences)
                if fig: st.plotly_chart(fig, use_container_width=True)
            with right:
                fig = plot_length_distribution(sequences)
                if fig: st.plotly_chart(fig, use_container_width=True)

            # Base composition if present
            if sequences and "A" in sequences[0] and "denom_used" in sequences[0]:
                fig = plot_base_composition(sequences)
                if fig: st.plotly_chart(fig, use_container_width=True)

            # Quality distribution if FASTQ fields exist
            qfig = plot_quality_distribution(sequences)
            if qfig is not None:
                st.plotly_chart(qfig, use_container_width=True)

        # ---- Tab 3: Details (sequence selector + truncated sequence) ----
        with tab3:
            st.subheader("Detailed Sequence Information")

            # Choose label field
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
                st.markdown(f"**GC%:** {gc_val if gc_val is not None else 'N/A'}")

                if "Avg_quality" in selected:
                    st.markdown(f"**Avg Quality:** {selected['Avg_quality']:.2f}")

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
                label="⬇️ Download as CSV",
                data=csv,
                file_name=f"{uploaded_file.name}_processed.csv",
                mime="text/csv",
            )
