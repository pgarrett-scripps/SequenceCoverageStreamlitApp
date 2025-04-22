from functools import lru_cache
import pandas as pd
import streamlit as st
import peptacular as pt

@lru_cache(maxsize=32)
def extract_peptides_from_file_cached(file_content, file_name):
    """Extract peptide sequences from file contents."""

    if file_name.endswith('.txt'):
        try:
            # Process directly from content without writing to disk
            lines = file_content.decode('utf-8').splitlines()
            peptides = [line.strip() for line in lines if line.strip()]
            return peptides
        except Exception as e:
            st.error(f"Error reading file {file_name}: {str(e)}")
            return []

    if file_name.endswith('.csv'):  # extract Sequence column
        try:
            # Process directly from content without writing to disk
            from io import BytesIO
            df = pd.read_csv(BytesIO(file_content))
            if 'Sequence' in df.columns:
                peptides = df['Sequence'].dropna().tolist()
                return peptides
            else:
                st.error("CSV file must contain a 'Sequence' column.")
                return []
        except Exception as e:
            st.error(f"Error reading file {file_name}: {str(e)}")
            return []


def extract_peptides_from_file(file):
    """Extract peptide sequences from an uploaded file."""
    return extract_peptides_from_file_cached(file.getvalue(), file.name)

@lru_cache(maxsize=32)
def calculate_coverage(protein_seq, peptides_tuple, accumulate):
    """Calculate coverage with caching for performance."""
    return pt.coverage(protein_seq, peptides_tuple, accumulate)