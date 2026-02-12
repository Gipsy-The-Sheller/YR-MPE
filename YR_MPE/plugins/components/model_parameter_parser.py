"""
Model Parameter Parser for IQ-TREE .iqtree files
Extracts Q matrix, state frequencies, gamma parameters, and invariant sites
"""

import re
from typing import Dict, List, Optional, Tuple
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def validate_dna_sequences(sequences: List[SeqRecord]) -> Tuple[bool, str]:
    """
    Validate that sequences are DNA sequences (not protein)
    
    Args:
        sequences: List of SeqRecord objects
        
    Returns:
        (is_valid, error_message)
    """
    if not sequences:
        return False, "No sequences provided"
    
    # Define valid DNA characters (including IUPAC degenerate bases)
    valid_dna_chars = set('ACGTURYKMSWBDHVNacgturykmswbdhvn-?!.')
    
    # Define protein-specific characters
    protein_chars = set('EFILPQZ*efilpqz*')
    
    for seq_record in sequences:
        seq_str = str(seq_record.seq)
        for char in seq_str:
            char_upper = char.upper()
            # Check for protein-specific characters
            if char_upper in protein_chars:
                return False, f"Protein-specific character '{char}' found in sequence '{seq_record.id}'. This plugin only supports DNA sequences."
            # Check for invalid characters
            if char_upper not in valid_dna_chars:
                return False, f"Invalid character '{char}' found in sequence '{seq_record.id}'. This plugin only supports DNA sequences (A, C, G, T and IUPAC degenerate bases)."
    
    return True, ""


def parse_iqtree_file(iqtree_file_path: str) -> Optional[Dict]:
    """
    Parse .iqtree file and extract model parameters
    
    Args:
        iqtree_file_path: Path to .iqtree file
        
    Returns:
        Dictionary containing model parameters or None if parsing fails
    """
    try:
        with open(iqtree_file_path, 'r', encoding='utf-8') as f:
            content = f.read()
        
        result = {
            'model_name': None,
            'rate_parameters': {},
            'state_frequencies': {},
            'q_matrix': [],
            'rate_heterogeneity': {
                'model': None,
                'gamma_alpha': None,
                'invariant_sites': None
            }
        }
        
        # Parse model name
        model_match = re.search(r'Model of substitution:\s*(.+)', content)
        if model_match:
            result['model_name'] = model_match.group(1).strip()
        
        # Parse rate parameters R
        rate_pattern = r'Rate parameter R:\s*\n\s*A-C:\s*([0-9.]+)\s*\n\s*A-G:\s*([0-9.]+)\s*\n\s*A-T:\s*([0-9.]+)\s*\n\s*C-G:\s*([0-9.]+)\s*\n\s*C-T:\s*([0-9.]+)\s*\n\s*G-T:\s*([0-9.]+)'
        rate_match = re.search(rate_pattern, content)
        if rate_match:
            result['rate_parameters'] = {
                'A-C': float(rate_match.group(1)),
                'A-G': float(rate_match.group(2)),
                'A-T': float(rate_match.group(3)),
                'C-G': float(rate_match.group(4)),
                'C-T': float(rate_match.group(5)),
                'G-T': float(rate_match.group(6))
            }
        
        # Parse state frequencies
        freq_pattern = r'State frequencies:.*?\n\s*pi\(A\)\s*=\s*([0-9.]+)\s*\n\s*pi\(C\)\s*=\s*([0-9.]+)\s*\n\s*pi\(G\)\s*=\s*([0-9.]+)\s*\n\s*pi\(T\)\s*=\s*([0-9.]+)'
        freq_match = re.search(freq_pattern, content, re.DOTALL)
        if freq_match:
            result['state_frequencies'] = {
                'pi(A)': float(freq_match.group(1)),
                'pi(C)': float(freq_match.group(2)),
                'pi(G)': float(freq_match.group(3)),
                'pi(T)': float(freq_match.group(4))
            }
        
        # Parse Q matrix
        q_pattern = r'Rate matrix Q:\s*\n\s*A\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s*\n\s*C\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s*\n\s*G\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s*\n\s*T\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)'
        q_match = re.search(q_pattern, content)
        if q_match:
            values = [float(q_match.group(i)) for i in range(1, 17)]
            result['q_matrix'] = [
                values[0:4],   # A row
                values[4:8],   # C row
                values[8:12],  # G row
                values[12:16]  # T row
            ]
        
        # Parse rate heterogeneity model
        hetero_match = re.search(r'Model of rate heterogeneity:\s*(.+)', content)
        if hetero_match:
            result['rate_heterogeneity']['model'] = hetero_match.group(1).strip()
            
            # Parse gamma shape alpha
            gamma_match = re.search(r'Gamma shape alpha:\s*([0-9.]+)', content)
            if gamma_match:
                result['rate_heterogeneity']['gamma_alpha'] = float(gamma_match.group(1))
            
            # Parse invariant sites proportion
            invariant_match = re.search(r'Proportion of invariable sites:\s*([0-9.]+)', content)
            if invariant_match:
                result['rate_heterogeneity']['invariant_sites'] = float(invariant_match.group(1))
        
        return result
        
    except Exception as e:
        print(f"Error parsing .iqtree file: {e}")
        return None


def get_q_matrix_color_scale(q_matrix: List[List[float]]) -> Tuple[float, float]:
    """
    Get the min and max values from Q matrix (excluding diagonal)
    for color scaling
    
    Args:
        q_matrix: 4x4 Q matrix
        
    Returns:
        (min_value, max_value)
    """
    all_values = []
    for i in range(4):
        for j in range(4):
            if i != j:  # Exclude diagonal
                all_values.append(q_matrix[i][j])
    
    if not all_values:
        return 0.0, 1.0
    
    return min(all_values), max(all_values)


def format_q_matrix_for_display(q_matrix: List[List[float]], precision: int = 4) -> List[List[str]]:
    """
    Format Q matrix values for display
    
    Args:
        q_matrix: 4x4 Q matrix
        precision: Number of decimal places
        
    Returns:
        Formatted matrix as list of lists of strings
    """
    formatted = []
    for row in q_matrix:
        formatted_row = [f"{value:.{precision}f}" for value in row]
        formatted.append(formatted_row)
    return formatted