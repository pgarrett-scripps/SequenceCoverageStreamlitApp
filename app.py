import streamlit as st
import pandas as pd
import numpy as np
import plotly.graph_objects as go
import peptacular as pt
from functools import lru_cache

st.set_page_config(layout="wide")


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

# App settings
st.sidebar.header("Settings")

# Visualization settings section
st.sidebar.subheader("Visualization Settings")
show_amino_acids = st.sidebar.toggle("Show amino acids on x-axis", value=False,
                                     help="Toggle between position numbers and amino acid letters")
show_data_editor = st.sidebar.toggle("Show data table (slow for large sequences)", value=False,
                                     help="Display raw data in table format")
show_file_labels = st.sidebar.toggle("Show file name labels", value=True,
                                     help="Show or hide file names above each coverage bar")
show_y_axis_labels = st.sidebar.toggle("Show y-axis labels", value=False,
                                       help="Show or hide file names on the y-axis")

label_font_size = st.sidebar.number_input("Label font size",
                                          min_value=8,
                                          max_value=50,
                                          value=12,
                                          step=1,
                                          help="Font size for file name labels")

# For long sequences, downsample the amino acid ticks
downsample_factor = st.sidebar.slider("X-axis label frequency",
                                      min_value=1,
                                      max_value=100,
                                      value=10,
                                      help="For long sequences, show every Nth position")

y_spacing = st.sidebar.number_input("Vertical spacing between bars",
                                    min_value=0.1,
                                    max_value=10.0,
                                    value=1.5,
                                    step=0.1,
                                    help="Adjust the vertical spacing between coverage bars")

label_color = st.sidebar.color_picker("Label color",
                                      value="#000000",
                                      help="Color for file name labels")

bg_color_rgba = st.sidebar.color_picker("Label background color",
                                        value="#FFFFFF",
                                        help="Background color for file name labels")

# Add colormap selection
colormap_options = [
    'Viridis', 'Plasma', 'Inferno', 'Magma', 'Cividis',
    'Turbo', 'Blues', 'Greens', 'Reds', 'Purples',
    'YlOrRd', 'YlGnBu', 'RdBu', 'Jet', 'Rainbow'
]
selected_colormap = st.sidebar.selectbox(
    "Color scheme",
    options=colormap_options,
    index=0,
    help="Select color scheme for coverage visualization"
)

st.caption('Motif')

motif_regex = st.sidebar.text_input(
    "Motif regex",
    value="N[ST][^P]",
    help="Regular expression to match specific motifs in the protein sequence."
)

motif_line_width = st.sidebar.number_input(
    "Motif line width",
    min_value=1,
    max_value=10,
    value=2,
    step=1,
    help="Width of the motif line in the visualization."
)
motif_color = st.sidebar.color_picker(
    "Motif line color",
    value="#FF0000",
    help="Color for the motif line in the visualization."
)

motif_text_size = st.sidebar.number_input(
    "Motif text size",
    min_value=8,
    max_value=50,
    value=12,
    step=1,
    help="Font size for the motif text in the visualization."
)


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

# Cached coverage calculation


@lru_cache(maxsize=32)
def calculate_coverage(protein_seq, peptides_tuple, accumulate):
    """Calculate coverage with caching for performance."""
    return pt.coverage(protein_seq, peptides_tuple, accumulate)


if not st.button("Analyze Coverage"):
    if not protein_sequence:
        st.error("Please enter a protein sequence.")
    elif not uploaded_files:
        st.error("Please upload at least one file with peptide sequences.")
    else:

        # Process each file
        coverage_data = {}
        coverage_percentages = {}
        for i, uploaded_file in enumerate(uploaded_files):
            # First half of progress
            peptides = extract_peptides_from_file(uploaded_file)

            if peptides:

                # Convert list to tuple for caching
                peptides_tuple = tuple(peptides)
                coverage_array = calculate_coverage(
                    protein_sequence, peptides_tuple, True)

                binary_coverage_array = calculate_coverage(
                    protein_sequence, peptides_tuple, False)

                if coverage_array is not None:
                    coverage_data[uploaded_file.name] = coverage_array
                else:
                    st.warning(
                        f"Could not calculate coverage for {uploaded_file.name}")

                if binary_coverage_array is not None:
                    coverage_percentages[uploaded_file.name] = np.sum(
                        binary_coverage_array) / len(binary_coverage_array) * 100
                else:
                    st.warning(
                        f"Could not calculate coverage for {uploaded_file.name}")

            else:
                st.warning(
                    f"No valid peptide sequences found in {uploaded_file.name}")

        motif_sites = []
        if motif_regex:
            # Find motif sites in the protein sequence
            motif_sites = list(pt.get_cleavage_sites(
                protein_sequence, motif_regex))

        motif_ranges = []
        if motif_regex:
            # Find motif ranges in the protein sequence
            motif_ranges = list(pt.get_regex_match_range(
                protein_sequence, motif_regex, 0))

        # allow users to updaet filenames:
        name_map = {}
        for fn in coverage_data.keys():
            new_name = st.text_input(
                f"Update filename for {fn}:",
                value=fn,
                help="Update the filename for better clarity."
            )
            name_map[fn] = new_name

        # Update coverage data with new names
        coverage_data = {name_map[fn]: coverage_data[fn]
                         for fn in coverage_data.keys()}

        # Update coverage percentages with new names
        coverage_percentages = {name_map[fn]: coverage_percentages[fn]
                                for fn in coverage_percentages.keys()}

        # Create a plotly figure for coverage visualization
        fig = go.Figure()

        # Get unique files for plotting
        files = sorted(coverage_data.keys())

        # Only create DataFrame if we're going to show it

        min_index, max_index = st.slider(
            "Select range of positions to view:",
            min_value=1,
            max_value=len(protein_sequence)-1,
            value=(1, len(protein_sequence)-1),
            step=1,
            help="Select the range of positions to view in the coverage data."
        )

        # remove oob values
        for i, file in enumerate(files):
            coverage_array = coverage_data[file]
            coverage_data[file] = coverage_array[min_index-1:max_index]

        # filter motif sites:
        motif_sites = [
            site for site in motif_sites if min_index <= site <= max_index]

        motif_ranges = [
            (start, end) for start, end in motif_ranges if min_index <= start <= max_index and min_index <= end <= max_index]

        cov_data = []
        for filename, coverage_array in coverage_data.items():
            for i, covered in enumerate(coverage_array, min_index-1):
                if covered > 0:  # Only include positions with coverage
                    cov_data.append({
                        'File': filename,
                        'index': i+1,
                        'amino_acid': protein_sequence[i],
                        'Coverage': covered
                    })

        if not cov_data:
            st.warning("No coverage data available for the selected range.")
            st.stop()

        # Create a DataFrame for coverage data
        cov_df = pd.DataFrame(cov_data)

        max_cov = cov_df['Coverage'].max()
        min_cov = cov_df['Coverage'].min()

        st.dataframe(cov_df, use_container_width=True)

        # Define y-positions for each file with adequate spacing
        y_positions = list(range(len(files)))
        y_positions.reverse()  # Reverse to keep same order as before

        # Add faint horizontal lines to indicate protein sequence length
        for i, file in enumerate(files):
            y_pos = y_positions[i] * y_spacing  # Use spacing to position bars

            # Add horizontal line for the protein length
            fig.add_shape(
                type="line",
                x0=min_index-1,  # Start at position 1
                x1=max_index+1,  # End at the length of the protein
                y0=y_pos,
                y1=y_pos,
                line=dict(
                    color="rgba(180, 180, 180, 0.3)",
                    width=2,
                    dash="solid"
                ),
                layer="below"  # Ensure the line is drawn behind everything
            )

        # add verticle lines for motif sites
        for site in motif_sites:
            # Add vertical line
            fig.add_shape(
                type="line",
                x0=site,
                x1=site,
                # Extend slightly above the top bar
                y0=-0.5,
                y1=(max(y_positions) + 0.7) * y_spacing,
                line=dict(
                    color=motif_color,
                    width=motif_line_width,
                    dash="solid"
                ),
                layer="above"  # Ensure the line is drawn above everything
            )

            # Add text annotation above the line with vertical text
            fig.add_annotation(
                x=site,
                # Position text above the top of the line
                y=(max(y_positions) + 0.8) * y_spacing,
                text=f"{site}",
                showarrow=False,
                font=dict(
                    size=motif_text_size,
                    color=motif_color,
                ),
                xanchor='center',
                yanchor='bottom',
                textangle=270  # Makes the text vertical (rotated 270 degrees)
            )

        # For each file, create a horizontal bar representation
        for i, file in enumerate(files):
            coverage_array = coverage_data[file]
            y_pos = y_positions[i] * y_spacing  # Use spacing to position bars

            # Group consecutive positions with the SAME coverage value - optimized approach
            segments = []
            start_pos = None
            end_pos = None
            current_coverage = None

            # Iterate through coverage array directly
            for pos, coverage in enumerate(coverage_array, min_index-1):
                if coverage > 0:
                    if start_pos is None:
                        # Start a new segment
                        start_pos = pos + 1  # 1-indexed
                        end_pos = pos + 1
                        current_coverage = coverage
                    elif pos == end_pos and coverage == current_coverage:
                        # Continue segment only if coverage value is the same
                        end_pos = pos + 1
                    else:
                        # Save completed segment and start new one
                        segments.append({
                            'start': start_pos,
                            'end': end_pos,
                            'coverage': current_coverage
                        })
                        start_pos = pos + 1
                        end_pos = pos + 1
                        current_coverage = coverage

            # Add the last segment if exists
            if start_pos is not None:
                segments.append({
                    'start': start_pos,
                    'end': end_pos,
                    'coverage': current_coverage
                })

            # Process segments for plotting
            coverage_values = []
            bar_widths = []
            bar_positions = []
            hover_texts = []

            for segment in segments:
                start_pos = segment['start']
                end_pos = segment['end']
                coverage_val = segment['coverage']

                bar_widths.append(end_pos - start_pos + 1)
                bar_positions.append(start_pos - 0.5)
                coverage_values.append(coverage_val)
                hover_texts.append(
                    f'File: {file}<br>Positions: {start_pos}-{end_pos}<br>Coverage: {coverage_val}')

            # Add a single trace for all segments in this file with numerical coloring
            fig.add_trace(go.Bar(
                x=bar_widths,
                y=[y_pos] * len(bar_widths),  # Use custom y position
                orientation='h',
                base=bar_positions,
                width=0.5,  # Added this line to make bars thinner
                marker=dict(
                    color=coverage_values,
                    colorscale=selected_colormap,  # Use the selected colormap
                    # Only show color scale for the last file
                    showscale=True if i == len(files) - 1 else False,
                    cmin=0,  # Set minimum color value
                    cmax=max_cov,  # Set maximum color value from global max
                    # Add transparency (0 is fully transparent, 1 is solid)
                    opacity=0.85,
                    colorbar=dict(
                        title=dict(
                            text="Coverage",
                            side="right"  # Make the title vertical
                        ),
                        thickness=15,
                        # Adjust length to match the number of files
                        len=1.2,
                        tickmode='array',  # Use array mode for ticks
                        # Integer ticks from 0 to max
                        tickvals=list(range(0, int(max_cov) + 1)),
                        ticktext=[str(i) for i in range(0, int(max_cov) + 1)]
                    ) if i == len(files) - 1 else None,
                ),
                text=None,  # Remove text from bars
                textposition='none',  # Ensure no text is shown
                hovertext=hover_texts,
                hoverinfo='text',
                showlegend=False,
            ))

            # Add centered filename above the bars - only if show_file_labels is True
            if show_file_labels:
                # Find the middle position of the protein sequence
                mid_pos = (max_index - min_index) / 2

                # Add annotation for the filename with customized appearance
                fig.add_annotation(
                    x=mid_pos,
                    y=y_pos,  # Position slightly above the bar
                    text=file,
                    showarrow=False,
                    font=dict(
                        size=label_font_size,
                        color=label_color,
                        weight='bold'  # Make the text bold
                    ),
                    xanchor='center',
                    bgcolor=bg_color_rgba  # Use the RGBA background color
                )

        # Create optimized ticktext and tickvals based on toggle selection
        if show_amino_acids:
            # For long sequences, downsample the ticks
            seq_len = max_index - min_index + 1
            # Only show every Nth position
            tickvals = list(range(min_index, max_index + 1, downsample_factor))
            ticktext = [protein_sequence[i-1] for i in tickvals]
        else:
            # Show position numbers
            ticktext = list(range(min_index, max_index + 1, downsample_factor))
            tickvals = list(range(min_index, max_index + 1, downsample_factor))

        # Update layout
        fig.update_layout(
            title="Protein Sequence Coverage",
            xaxis_title="Amino Acid Position",
            height=max(300, 150 * len(files)),  # Adjusted for spacing
            showlegend=False,
            xaxis=dict(
                range=[min_index - 2, max_index + 2],
                ticktext=ticktext,
                tickvals=tickvals,
                showgrid=True,
                gridcolor='rgba(128, 128, 128, 0.5)',
                gridwidth=1,
                griddash='dot',
                autorange=True
            ),
            barmode='overlay',
            bargap=0,  # Remove default gap between bars
            yaxis=dict(
                showticklabels=show_y_axis_labels,  # Show y-axis labels based on toggle
                # Set tick mode to array when showing labels
                tickmode='array' if show_y_axis_labels else None,
                # Position of each file's bar
                tickvals=[
                    y_positions[i] * y_spacing for i in range(len(files))] if show_y_axis_labels else None,
                ticktext=files if show_y_axis_labels else None,  # Use file names as labels
                tickfont=dict(
                    size=label_font_size,
                    color=label_color
                ) if show_y_axis_labels else None,
                showgrid=False,
                # Add these parameters to move labels closer to the plot
            ),
        )

        # Create a secondary y-axis trace to show coverage percentages
        for i, file in enumerate(files):
            percentage = coverage_percentages[file]
            y_pos = y_positions[i] * y_spacing

            # Add a text annotation for the percentage on the right side
            fig.add_annotation(
                x=max_index + 2,  # Position slightly to the right of the maximum x value
                y=y_pos,
                text=f"{percentage:.1f}%",
                showarrow=False,
                font=dict(
                    size=label_font_size,
                    color=label_color,
                ),
                xanchor='left',
                yanchor='middle',
                bgcolor=bg_color_rgba  # Use the same background color as file labels
            )

        # After you've defined motif_ranges and before creating the figure:
        for start, end in motif_ranges:
            # Adjust coordinates to account for min_index filtering
            adj_start = max(start, min_index)
            adj_end = min(end, max_index)

            # Only add shapes for motifs that fall within the visible range
            if adj_start <= max_index and adj_end >= min_index:
                fig.add_shape(
                    type="rect",
                    x0=adj_start + 1 - 0.5,  # Offset by 0.5 to align with bar positions
                    x1=adj_end + 0.5,
                    y0=-0.5,  # Extend slightly below the bottom
                    # Match the motif line height
                    y1=(max(y_positions) + 0.7) * y_spacing,
                    # Light color for background
                    fillcolor="rgba(255, 235, 205, 0.85)",
                    line=dict(width=0),  # No border
                    layer="below"  # Draw below all other elements
                )

        # Display the figure
        st.subheader("Coverage Visualization")
        st.plotly_chart(fig, use_container_width=True)
