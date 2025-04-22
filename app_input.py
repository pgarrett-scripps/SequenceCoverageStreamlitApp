import numpy as np
import pandas as pd
import streamlit as st
from dataclasses import dataclass
from typing import List, Optional, Tuple, Any

import peptacular as pt

from util import calculate_coverage, extract_peptides_from_file

DEFAULT_PROTIEN = """MFVFLVLLPLVSSQCVNLTTRTQLPPAYTNSFTRGVYYPDKVFRSSVLHSTQDLFLPFFSNVTWFHAIHVSGTNGTKRFDNPVLPFNDGVYFASTEKSNIIRGWIFGTTLDSKTQSLLIVNNATNVVIKVCEFQFCNDPFLGVYYHKNNKSWMESEFRVYSSANNCTFEYVSQPFLMDLEGKQGNFKNLREFVFKNIDGYFKIYSKHTPINLVRDLPQGFSALEPLVDLPIGINITRFQTLLALHRSYLTPGDSSSGWTAGAAAYYVGYLQPRTFLLKYNENGTITDAVDCALDPLSETKCTLKSFTVEKGIYQTSNFRVQPTESIVRFPNITNLCPFGEVFNATRFASVYAWNRKRISNCVADYSVLYNSASFSTFKCYGVSPTKLNDLCFTNVYADSFVIRGDEVRQIAPGQTGKIADYNYKLPDDFTGCVIAWNSNNLDSKVGGNYNYLYRLFRKSNLKPFERDISTEIYQAGSTPCNGVEGFNCYFPLQSYGFQPTNGVGYQPYRVVVLSFELLHAPATVCGPKKSTNLVKNKCVNFNFNGLTGTGVLTESNKKFLPFQQFGRDIADTTDAVRDPQTLEILDITPCSFGGVSVITPGTNTSNQVAVLYQDVNCTEVPVAIHADQLTPTWRVYSTGSNVFQTRAGCLIGAEHVNNSYECDIPIGAGICASYQTQTNSPRRARSVASQSIIAYTMSLGAENSVAYSNNSIAIPTNFTISVTTEILPVSMTKTSVDCTMYICGDSTECSNLLLQYGSFCTQLNRALTGIAVEQDKNTQEVFAQVKQIYKTPPIKDFGGFNFSQILPDPSKPSKRSFIEDLLFNKVTLADAGFIKQYGDCLGDIAARDLICAQKFNGLTVLPPLLTDEMIAQYTSALLAGTITSGWTFGAGAALQIPFAMQMAYRFNGIGVTQNVLYENQKLIANQFNSAIGKIQDSLSSTASALGKLQDVVNQNAQALNTLVKQLSSNFGAISSVLNDILSRLDKVEAEVQIDRLITGRLQSLQTYVTQQLIRAAEIRASANLAATKMSECVLGQSKRVDFCGKGYHLMSFPQSAPHGVVFLHVTYVPAQEKNFTTAPAICHDGKAHFPREGVFVSNGTHWFVTQRNFYEPQIITTDNTFVSGNCDVVIGIVNNTVYDPLQPELDSFKEELDKYFKNHTSPDVDLGDISGINASVVNIQKEIDRLNEVAKNLNESLIDLQELGKYEQYIKWPWYIWLGFIAGLIAIVMVTIMLCCMTSCCSCLKGCCSCGSCCKFDEDDSEPVLKGVKLHYT"""


@dataclass
class AppInputs:
    # Sequence data
    protein_sequence: str
    uploaded_files: List[Any]

    # Visualization settings
    show_amino_acids: bool
    show_file_labels: bool
    show_y_axis_labels: bool
    label_font_size: int
    downsample_factor: int
    y_spacing: float
    label_color: str
    bg_color_rgba: str
    selected_colormap: str
    min_coverage: int
    binary_coverage: bool
    display_zero_segments: bool

    # Motif settings
    motif_regex: str
    motif_line_width: int
    motif_color: str
    motif_text_size: int
    motif_background_color: str
    show_motif_ranges: bool
    show_motif_sites: bool

    min_index: int
    max_index: int
    name_map: dict
    custom_sites: List[int]

    @property
    def motif_sites(self) -> List[int]:
        motif_sites = []
        if self.motif_regex:
            motif_sites = list(pt.get_cleavage_sites(
                self.protein_sequence, self.motif_regex))

        # Filter motif sites based on min_index and max_index
        motif_sites = [
            site for site in motif_sites if self.min_index <= site <= self.max_index]

        return motif_sites

    @property
    def motif_ranges(self) -> List[Tuple[int, int]]:

        motif_ranges = []
        if self.motif_regex:
            motif_ranges = list(pt.get_regex_match_range(
                self.protein_sequence, self.motif_regex, 0))

        # Filter motif ranges based on min_index and max_index
        motif_ranges = [
            (start, end) for start, end in motif_ranges if self.min_index <= start <= self.max_index]

        return motif_ranges

    @property
    def coverage_data(self) -> dict:
        # Process each file
        coverage_data = {}
        for i, uploaded_file in enumerate(self.uploaded_files):
            # First half of progress
            peptides = extract_peptides_from_file(uploaded_file)

            if peptides:

                # Convert list to tuple for caching
                peptides_tuple = tuple(peptides)
                coverage_array = calculate_coverage(
                    self.protein_sequence, peptides_tuple, True)

                if coverage_array is not None:
                    coverage_data[self.name_map[uploaded_file.name]
                                  ] = coverage_array
                else:
                    st.warning(
                        f"Could not calculate coverage for {uploaded_file.name}")

            else:
                st.warning(
                    f"No valid peptide sequences found in {uploaded_file.name}")

        # filter min and max
        for file_name, coverage_array in coverage_data.items():
            # Filter the coverage array based on min_index and max_index
            filtered_coverage = coverage_array[self.min_index-1:self.max_index]
            coverage_data[file_name] = filtered_coverage

        # filter based on min_coverage
        for file_name, coverage_array in coverage_data.items():
            # Convert to numpy array and filter based on min_coverage
            coverage_array = np.array(coverage_array)
            filtered_coverage = np.where(
                coverage_array >= self.min_coverage, coverage_array, 0)
            coverage_data[file_name] = filtered_coverage

        # convert to binary if binary_coverage is selected
        if self.binary_coverage:
            for file_name, coverage_array in coverage_data.items():
                # Convert the coverage array to binary (0 or 1)
                binary_coverage = np.where(coverage_array > 0, 1, 0)
                coverage_data[file_name] = binary_coverage

        return coverage_data

    @property
    def coverage_percentages_dict(self) -> dict:
        coverage_date = self.coverage_data
        coverage_percentages = {}
        for file_name, coverage_array in coverage_date.items():
            # Calculate the percentage coverage
            total_positions = len(coverage_array)
            covered_positions = np.count_nonzero(coverage_array)
            coverage_percentage = (covered_positions / total_positions) * 100
            coverage_percentages[file_name] = coverage_percentage

        return coverage_percentages

    @property
    def cov_df(self) -> pd.DataFrame:
        cov_data = []
        for filename, coverage_array in self.coverage_data.items():
            for i, covered in enumerate(coverage_array, self.min_index-1):
                if covered > 0:  # Only include positions with coverage
                    cov_data.append({
                        'File': filename,
                        'index': i+1,
                        'amino_acid': self.protein_sequence[i],
                        'Coverage': covered
                    })

        # Create a DataFrame from the coverage data
        coverage_df = pd.DataFrame(cov_data)

        if len(coverage_df) == 0:
            return pd.DataFrame({'File': [], 'index': [], 'amino_acid': [], 'Coverage': []})

        return coverage_df


def setup_inputs() -> AppInputs:
    """Set up all inputs and return as a dataclass."""

    # Input for protein sequence
    st.header("Protein Sequence")
    protein_sequence = st.text_area(
        "Enter protein sequence:",
        value=DEFAULT_PROTIEN,
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
    st.header("Settings")

    # setup tabs
    tabs = st.tabs(["Format", "X-Axis", "Motif", "Names"])

    with tabs[0]:
        show_file_labels = st.toggle("Overlay file lables", value=False,
                                     help="Show or hide file names above each coverage bar")
        show_y_axis_labels = st.toggle("Show Y-axis files", value=True,
                                       help="Show or hide file names on the y-axis")

        label_font_size = st.number_input("Label font size",
                                          min_value=8,
                                          max_value=50,
                                          value=12,
                                          step=1,
                                          help="Font size for file name labels")

        y_spacing = st.number_input("Vertical spacing between bars",
                                    min_value=0.1,
                                    max_value=20.0,
                                    value=1.5,
                                    step=0.1,
                                    help="Adjust the vertical spacing between coverage bars")

        c1, c2 = st.columns(2)
        label_color = c1.color_picker("Label color",
                                      value="#000000",
                                      help="Color for file name labels")
        bg_color_rgba = c2.color_picker("Label background color",
                                        value="#FFFFFF",
                                        help="Background color for file name labels")

        # Add colormap selection
        colormap_options = ['aggrnyl', 'agsunset', 'algae', 'amp', 'armyrose', 'balance', 'blackbody', 'bluered', 'blues', 'blugrn', 'bluyl', 'brbg', 'brwnyl', 'bugn', 'bupu', 'burg', 'burgyl', 'cividis', 'curl', 'darkmint', 'deep', 'delta', 'dense', 'earth', 'edge', 'electric', 'emrld', 'fall', 'geyser', 'gnbu', 'gray', 'greens', 'greys', 'haline', 'hot', 'hsv', 'ice', 'icefire', 'inferno', 'jet', 'magenta', 'magma', 'matter', 'mint', 'mrybm', 'mygbm', 'oranges',
                            'orrd', 'oryel', 'oxy', 'peach', 'phase', 'picnic', 'pinkyl', 'piyg', 'plasma', 'plotly3', 'portland', 'prgn', 'pubu', 'pubugn', 'puor', 'purd', 'purp', 'purples', 'purpor', 'rainbow', 'rdbu', 'rdgy', 'rdpu', 'rdylbu', 'rdylgn', 'redor', 'reds', 'solar', 'spectral', 'speed', 'sunset', 'sunsetdark', 'teal', 'tealgrn', 'tealrose', 'tempo', 'temps', 'thermal', 'tropic', 'turbid', 'turbo', 'twilight', 'viridis', 'ylgn', 'ylgnbu', 'ylorbr', 'ylorrd']
        selected_colormap = st.selectbox(
            "Color scheme",
            options=colormap_options,
            index=colormap_options.index('cividis'),
            help="Select color scheme for coverage visualization"
        )

        min_coverage = st.number_input(
            "Minimum coverage",
            min_value=0,
            value=1,
            step=1,
            help="Minimum coverage value to display in the visualization."
        )

        binary_coverage = st.toggle("Binary coverage", value=False,
                                    help="Display coverage as binary (0 or 1) instead of percentage")

        display_zero_segments = st.toggle("Display zero segments", value=False,
                                          help="Display segments with zero coverage in the visualization")

    with tabs[1]:

        min_index, max_index = st.slider(
            "Select range of positions to view:",
            min_value=1,
            max_value=len(protein_sequence),
            value=(1, len(protein_sequence)),
            step=1,
            help="Select the range of positions to view in the coverage data."
        )

        # Visualization settings section
        show_amino_acids = st.toggle("Show amino acids on x-axis", value=False,
                                     help="Toggle between position numbers and amino acid letters")

        # For long sequences, downsample the amino acid ticks
        downsample_factor = st.number_input(
            "Downsample factor",
            min_value=1,
            max_value=len(protein_sequence),
            value=20,
            step=1,
            help="Downsample factor for amino acid ticks on the x-axis."
        )

    with tabs[2]:

        c1, c2 = st.columns([3, 1])
        motif_regex = c1.text_input(
            "Motif regex",
            value="N[ST][^P]",
            help="Regular expression to match specific motifs in the protein sequence."
        )

        motif_color = c2.color_picker(
            "Motif color",
            value="#FF0000",
            help="Color for the motif line in the visualization."
        )

        c1, c2 = st.columns(2)
        motif_line_width = c1.number_input(
            "Motif line width",
            min_value=1,
            max_value=10,
            value=2,
            step=1,
            help="Width of the motif line in the visualization."
        )

        motif_text_size = c2.number_input(
            "Motif text size",
            min_value=8,
            max_value=50,
            value=12,
            step=1,
            help="Font size for the motif text in the visualization."
        )

        motif_background_color = st.color_picker(
            "Motif background color",
            value="#FF0000",  # red
            help="Background color for the motif text in the visualization."
        )

        show_motif_ranges = st.toggle("Show motif ranges", value=False,
                                      help="Show or hide the ranges of the motifs in the visualization")

        show_motif_sites = st.toggle("Show motif sites", value=False,
                                     help="Show or hide the sites of the motifs in the visualization")

        custom_sites = st.text_input(
            "Custom sites",
            value="",
            help="Comma-separated list of custom sites to highlight in the visualization."
        )

        if custom_sites:
            try:
                custom_sites = [int(site.strip())
                                for site in custom_sites.split(",")]
                # Ensure the sites are within the range of the protein sequence
                custom_sites = [
                    site for site in custom_sites if 1 <= site <= len(protein_sequence)]
            except ValueError:
                st.error(
                    "Invalid input. Please enter a comma-separated list of integers.")
                custom_sites = []
        else:
            custom_sites = []

    with tabs[3]:
        name_map = {}
        # allow user to change file names
        st.subheader("File Names")
        st.write("Change the names of the uploaded files:")
        for i, uploaded_file in enumerate(uploaded_files):
            new_name = st.text_input(
                f"File {i+1} name",
                value=uploaded_file.name,
                key=f"file_name_{i}"
            )
            # Update the file name in the uploaded_files list
            name_map[uploaded_file.name] = new_name

    # Return all inputs as a dataclass
    return AppInputs(
        protein_sequence=protein_sequence,
        uploaded_files=uploaded_files,
        show_amino_acids=show_amino_acids,
        show_file_labels=show_file_labels,
        show_y_axis_labels=show_y_axis_labels,
        label_font_size=label_font_size,
        downsample_factor=downsample_factor,
        y_spacing=y_spacing,
        label_color=label_color,
        bg_color_rgba=bg_color_rgba,
        selected_colormap=selected_colormap,
        motif_regex=motif_regex,
        motif_line_width=motif_line_width,
        motif_color=motif_color,
        motif_text_size=motif_text_size,
        motif_background_color=motif_background_color,
        show_motif_ranges=show_motif_ranges,
        show_motif_sites=show_motif_sites,
        min_index=min_index,
        max_index=max_index,
        name_map=name_map,
        min_coverage=min_coverage,
        binary_coverage=binary_coverage,
        display_zero_segments=display_zero_segments,
        custom_sites=custom_sites
    )


# This line allows importing the setup_inputs function and AppInputs class
# from this module while also allowing the file to be run directly
if __name__ == "__main__":
    inputs = setup_inputs()
