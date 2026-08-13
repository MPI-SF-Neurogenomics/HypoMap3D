
import gzip
import re

class File_Analysis:
    """
    A class for analyzing sequences in a FASTQ file to determine the longest sequence length and T-richness (t_score).
    
    Attributes:
        filename (str): The name of the FASTQ file being analyzed.
        line_length (int): The length of the longest sequence in the file.
        t_score (float): A score representing the T-richness of sequences in the file.
    """
    def __init__(self, filename: str, max_lines: int = 36) -> None:
        """
        Initializes the file_analysis class by reading sequences from the given FASTQ file and calculating
        the line length and t_score.

        Args:
            filename (str): The name of the FASTQ file to analyze.
            max_lines (int, optional): The maximum number of sequences to analyze. Defaults to 36.
        """
        self.filename = filename
        
        sequences = []
        with gzip.open(filename, 'rt') as fp:
            for i, line in enumerate(fp):
                if i % 4 == 1:
                    sequences.append(line.strip())
                if i > max_lines:
                    break
        
        self.line_length = max(len(seq) for seq in sequences)
        
        total_sequence = "".join(sequences)
        t_matches = re.findall(r'T+', total_sequence.upper())
        
        if t_matches:
            self.t_score = 1 / sum(len(match) ** 2 for match in t_matches)
        else:
            self.t_score = 1

    def __repr__(self) -> str:
        """
        Returns a string representation of the file_analysis object, including the filename, t_score, and line_length.

        Returns:
            str: A string representation of the file_analysis object.
        """
        return f"{self.filename}\n\tt_score = {self.t_score}\n\tline_length = {self.line_length}"