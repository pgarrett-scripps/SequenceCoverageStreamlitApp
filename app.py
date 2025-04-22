import streamlit as st
import plotly.graph_objects as go

from app_input import setup_inputs

st.set_page_config(layout="wide")


st.title("Peptide Sequence Coverage Analyzer")

st.markdown(
    """
    <style>
        section[data-testid="stSidebar"] {
            width: 600px !important; # Set the width to your desired value
        }
    </style>
    """,
    unsafe_allow_html=True,
)

with st.sidebar:
    params = setup_inputs()


if not params.protein_sequence:
    st.error("Please enter a protein sequence.")
elif not params.uploaded_files:
    st.error("Please upload at least one file with peptide sequences.")
else:

    # Create a plotly figure for coverage visualization
    fig = go.Figure()

    # Get unique files for plotting
    files = sorted(params.coverage_data.keys())

    max_cov = params.cov_df['Coverage'].max()
    min_cov = params.cov_df['Coverage'].min()

    # Define y-positions for each file with adequate spacing
    y_positions = list(range(len(files)))
    y_positions.reverse()  # Reverse to keep same order as before

    # Add faint horizontal lines to indicate protein sequence length
    for i, file in enumerate(files):
        # Use spacing to position bars
        y_pos = y_positions[i] * params.y_spacing

        # Add horizontal line for the protein length
        fig.add_shape(
            type="line",
            x0=params.min_index-1,  # Start at position 1
            x1=params.max_index+1,  # End at the length of the protein
            y0=y_pos,
            y1=y_pos,
            line=dict(
                color="rgba(180, 180, 180, 0.3)",
                width=2,
                dash="solid"
            ),
            layer="below"  # Ensure the line is drawn behind everything
        )

    if params.show_motif_sites:
        # add verticle lines for motif sites
        for site in params.motif_sites:
            # Add vertical line
            fig.add_shape(
                type="line",
                x0=site,
                x1=site,
                # Extend slightly above the top bar
                y0=-0.5,
                y1=(max(y_positions) + 0.7) * params.y_spacing,
                line=dict(
                    color=params.motif_color,
                    width=params.motif_line_width,
                    dash="solid"
                ),
                layer="above"  # Ensure the line is drawn above everything
            )

            # Add text annotation above the line with vertical text
            fig.add_annotation(
                x=site,
                # Position text above the top of the line
                y=(max(y_positions) + 0.8) * params.y_spacing,
                text=f"{site}",
                showarrow=False,
                font=dict(
                    size=params.motif_text_size,
                    color=params.motif_color,
                ),
                xanchor='center',
                yanchor='bottom',
                textangle=270  # Makes the text vertical (rotated 270 degrees)
            )

    # For each file, create a horizontal bar representation
    for i, file in enumerate(files):
        coverage_array = params.coverage_data[file]
        # Use spacing to position bars
        y_pos = y_positions[i] * params.y_spacing

        # Group consecutive positions with the SAME coverage value - optimized approach
        segments = []
        start_pos = None
        end_pos = None
        current_coverage = None

        # Iterate through coverage array directly
        for pos, coverage in enumerate(coverage_array, params.min_index-1):
            # Only consider positions with coverage > 0
            if (params.display_zero_segments and coverage >= 0) or (params.display_zero_segments is False and coverage > 0):
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
                colorscale=params.selected_colormap,  # Use the selected colormap
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
        if params.show_file_labels:
            # Find the middle position of the protein sequence
            mid_pos = (params.max_index - params.min_index) / 2

            # Add annotation for the filename with customized appearance
            fig.add_annotation(
                x=mid_pos,
                y=y_pos,  # Position slightly above the bar
                text=file,
                showarrow=False,
                font=dict(
                    size=params.label_font_size,
                    color=params.label_color,
                    weight='bold'  # Make the text bold
                ),
                xanchor='center',
                bgcolor=params.bg_color_rgba  # Use the RGBA background color
            )

    # Create optimized ticktext and tickvals based on toggle selection
    if params.show_amino_acids:
        # For long sequences, downsample the ticks
        seq_len = params.max_index - params.min_index + 1
        # Only show every Nth position
        tickvals = list(
            range(params.min_index, params.max_index + 1, params.downsample_factor))
        ticktext = [params.protein_sequence[i-1] for i in tickvals]
    else:
        # Show position numbers
        ticktext = list(
            range(params.min_index, params.max_index + 1, params.downsample_factor))
        tickvals = list(
            range(params.min_index, params.max_index + 1, params.downsample_factor))

    # Update layout
    fig.update_layout(
        title="Protein Sequence Coverage",
        xaxis_title="Amino Acid Position",
        height=max(300, 150 * len(files)),  # Adjusted for spacing
        showlegend=False,
        xaxis=dict(
            range=[params.min_index - 2, params.max_index + 2],
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
            showticklabels=params.show_y_axis_labels,  # Show y-axis labels based on toggle
            # Set tick mode to array when showing labels
            tickmode='array' if params.show_y_axis_labels else None,
            # Position of each file's bar
            tickvals=[
                y_positions[i] * params.y_spacing for i in range(len(files))] if params.show_y_axis_labels else None,
            ticktext=files if params.show_y_axis_labels else None,  # Use file names as labels
            tickfont=dict(
                size=params.label_font_size,
                color=params.label_color
            ) if params.show_y_axis_labels else None,
            showgrid=False,
            # Add these parameters to move labels closer to the plot
        ),
    )

    # Create a secondary y-axis trace to show coverage percentages
    for i, file in enumerate(files):
        percentage = params.coverage_percentages_dict[file]
        y_pos = y_positions[i] * params.y_spacing

        # Add a text annotation for the percentage on the right side
        fig.add_annotation(
            x=params.max_index + 2,  # Position slightly to the right of the maximum x value
            y=y_pos,
            text=f"{percentage:.1f}%",
            showarrow=False,
            font=dict(
                size=params.label_font_size,
                color=params.label_color,
            ),
            xanchor='left',
            yanchor='middle',
            bgcolor=params.bg_color_rgba  # Use the same background color as file labels
        )

    # After you've defined motif_ranges and before creating the figure:
    if params.show_motif_ranges:
        for start, end in params.motif_ranges:
            # Adjust coordinates to account for min_index filtering
            adj_start = max(start, params.min_index)
            adj_end = min(end, params.max_index)

            # Only add shapes for motifs that fall within the visible range
            if adj_start <= params.max_index and adj_end >= params.min_index:
                fig.add_shape(
                    type="rect",
                    x0=adj_start + 1 - 0.5,  # Offset by 0.5 to align with bar positions
                    x1=adj_end + 0.5,
                    y0=-0.5,  # Extend slightly below the bottom
                    # Match the motif line height
                    y1=(max(y_positions) + 0.7) * params.y_spacing,
                    # Light color for background
                    fillcolor=params.motif_background_color,
                    # transparency
                    opacity=0.3,
                    line=dict(width=0),  # No border
                    layer="below"  # Draw below all other elements,
                )

    # show custom motif sites
    for site in params.custom_sites:
        # Add vertical line
        fig.add_shape(
            type="line",
            x0=site,
            x1=site,
            # Extend slightly above the top bar
            y0=-0.5,
            y1=(max(y_positions) + 0.7) * params.y_spacing,
            line=dict(
                color=params.motif_color,
                width=params.motif_line_width,
                dash="solid"
            ),
            layer="above"  # Ensure the line is drawn above everything
        )

        # Add text annotation above the line with vertical text
        fig.add_annotation(
            x=site,
            # Position text above the top of the line
            y=(max(y_positions) + 0.8) * params.y_spacing,
            text=f"{site}",
            showarrow=False,
            font=dict(
                size=params.motif_text_size,
                color=params.motif_color,
            ),
            xanchor='center',
            yanchor='bottom',
            textangle=270  # Makes the text vertical (rotated 270 degrees)
        )

    # Display the figure
    st.subheader("Coverage Visualization")
    st.plotly_chart(fig, use_container_width=True)
