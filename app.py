import streamlit as st
import pandas as pd
import numpy as np
import plotly.graph_objects as go
import peptacular.sequence.sequence_funcs as sequence_funcs
import tempfile
import os

st.title("Peptide Sequence Coverage Analyzer")

# Input for protein sequence
st.header("Protein Sequence")
protein_sequence = st.text_area(
    "Enter protein sequence:",
    height=100,
    help="Enter the protein amino acid sequence for coverage analysis."
)

# File uploader for peptide sequence files
st.header("Peptide Sequence Files")
uploaded_files = st.file_uploader(
    "Upload files containing peptide sequences (one per line):",
    accept_multiple_files=True,
    type=["txt", "csv", "fasta"],
    help="Each file should contain one peptide sequence per line."
)

def extract_peptides_from_file(file):
    """Extract peptide sequences from an uploaded file."""
    with tempfile.NamedTemporaryFile(delete=False) as temp_file:
        temp_file.write(file.getvalue())
        temp_file_path = temp_file.name
    
    try:
        with open(temp_file_path, 'r') as f:
            peptides = [line.strip() for line in f if line.strip()]
        os.unlink(temp_file_path)
        return peptides
    except Exception as e:
        os.unlink(temp_file_path)
        st.error(f"Error reading file {file.name}: {str(e)}")
        return []

def calculate_coverage_matrix(protein_seq, peptide_list):
    """
    Calculate the sequence coverage matrix for a list of peptides against a protein sequence.
    Returns a coverage array where 1 indicates covered position, 0 indicates not covered.
    """
    if not protein_seq or not peptide_list:
        return None
    
    coverage = np.zeros(len(protein_seq))
    
    for peptide in peptide_list:
        # Find all occurrences of the peptide in the protein sequence
        start_idx = 0
        while start_idx < len(protein_seq):
            idx = protein_seq.find(peptide, start_idx)
            if idx == -1:
                break
            # Mark these positions as covered
            coverage[idx:idx+len(peptide)] = 1
            start_idx = idx + 1
    
    return coverage

if st.button("Analyze Coverage"):
    if not protein_sequence:
        st.error("Please enter a protein sequence.")
    elif not uploaded_files:
        st.error("Please upload at least one file with peptide sequences.")
    else:
        st.subheader("Sequence Coverage Analysis")
        
        # Process each file
        coverage_data = {}
        for uploaded_file in uploaded_files:
            peptides = extract_peptides_from_file(uploaded_file)
            
            if peptides:
                st.write(f"File: {uploaded_file.name}, Number of peptides: {len(peptides)}")
                coverage_array = calculate_coverage_matrix(protein_sequence, peptides)
                if coverage_array is not None:
                    coverage_data[uploaded_file.name] = coverage_array
                    coverage_percentage = np.mean(coverage_array) * 100
                    st.write(f"Coverage: {coverage_percentage:.2f}%")
                else:
                    st.warning(f"Could not calculate coverage for {uploaded_file.name}")
        
        if coverage_data:
            # Create visualization
            fig = go.Figure()
            
            # X-axis will be amino acid positions
            x_positions = list(range(1, len(protein_sequence) + 1))
            
            # Add a trace for each file
            for i, (filename, coverage_array) in enumerate(coverage_data.items()):
                y_position = i
                
                # Create horizontal bars for covered regions
                segments = []
                start = None
                
                for pos, covered in enumerate(coverage_array):
                    if covered == 1 and start is None:
                        start = pos
                    elif covered == 0 and start is not None:
                        segments.append((start, pos-1))
                        start = None
                
                # Add final segment if needed
                if start is not None:
                    segments.append((start, len(coverage_array)-1))
                
                # Add bars for each segment
                for seg_start, seg_end in segments:
                    fig.add_trace(go.Bar(
                        x=[seg_end - seg_start + 1],
                        y=[filename],
                        orientation='h',
                        base=[seg_start],
                        name=filename,
                        showlegend=False,
                        marker_color=f'hsl({(i * 50) % 360}, 70%, 50%)'
                    ))
            
            # Configure layout
            fig.update_layout(
                title="Protein Sequence Coverage by Peptide File",
                xaxis_title="Amino Acid Position",
                yaxis_title="File",
                barmode='overlay',
                height=400 + (len(coverage_data) * 50),
                xaxis=dict(range=[0, len(protein_sequence)]),
                yaxis=dict(
                    categoryorder='array',
                    categoryarray=list(coverage_data.keys())
                )
            )
            
            st.plotly_chart(fig, use_container_width=True)
            
            # Show amino acid positions
            st.subheader("Protein Sequence with Position Labels")
            sequence_display = ""
            for i, aa in enumerate(protein_sequence):
                if i % 10 == 0:
                    sequence_display += f"<b>{i+1}</b>"
                sequence_display += aa
                if (i+1) % 50 == 0:
                    sequence_display += "<br>"
            
            st.markdown(f"<div style='font-family: monospace;'>{sequence_display}</div>", unsafe_allow_html=True)