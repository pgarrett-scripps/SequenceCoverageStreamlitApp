# Peptide Sequence Coverage Analyzer

A Streamlit application that analyzes peptide sequence coverage against a protein sequence and visualizes the results.

## Description

This application allows you to:

- Input a protein sequence
- Upload one or more files containing peptide sequences (one per line)
- Analyze the coverage of peptides against the protein sequence
- Visualize the coverage using an interactive Plotly chart

The application displays horizontal bars representing the covered regions of the protein for each uploaded peptide file, making it easy to visualize which parts of the protein are covered by different sets of peptides.

## Features

- Simple and intuitive web interface
- Support for multiple peptide sequence files
- Interactive visualization of sequence coverage
- Coverage percentage calculation for each file
- Protein sequence display with position labels

## Installation

### Prerequisites

- Python 3.7 or higher

### Setup

1. Clone this repository:
   ```
   git clone https://github.com/yourusername/SequenceCoverageStreamlitApp.git
   cd SequenceCoverageStreamlitApp
   ```

2. Install the required packages:
   ```
   pip install -r requirements.txt
   ```

## Usage

1. Run the Streamlit application:
   ```
   streamlit run app.py
   ```

2. Access the application in your web browser (typically at http://localhost:8501).

3. Enter your protein sequence in the text area.

4. Upload one or more files containing peptide sequences (one sequence per line).

5. Click "Analyze Coverage" to process the data and view the results.

## Example Peptide File Format

Each file should contain one peptide sequence per line, for example:

```
MLPDGR
YISFTQK
SPFLLR
AEYEPETLAK
```

## Technologies Used

- [Streamlit](https://streamlit.io/) - Web interface
- [Peptacular](https://peptacular.readthedocs.io/) - Peptide sequence analysis
- [Plotly](https://plotly.com/) - Interactive visualization
- [NumPy](https://numpy.org/) - Numerical operations
- [Pandas](https://pandas.pydata.org/) - Data manipulation

## License

This project is licensed under the MIT License - see the LICENSE file for details.