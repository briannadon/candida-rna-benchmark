#!/usr/bin/env python3
"""
Script to convert FASTA headers from ENA format to GFF3-compatible format.
ENA format: >ENA|CP017623|CP017623.1 Candida albicans SC5314 chromosome 1 sequence.
GFF3 format: >Ca22chr1A_C_albicans_SC5314
"""

import sys
import re

def convert_header(header):
    """
    Convert ENA FASTA header to GFF3-compatible format.
    
    Args:
        header (str): FASTA header line starting with '>'
        
    Returns:
        str: Converted header in GFF3 format
    """
    # Pattern to extract chromosome number/letter
    pattern = r'chromosome\s+([0-9R]+)'
    match = re.search(pattern, header)
    
    if not match:
        # If no chromosome found, return original header
        return header
    
    chrom = match.group(1)
    
    # Handle special case for chromosome R
    if chrom == 'R':
        return ">Ca22chrRA_C_albicans_SC5314"
    else:
        return f">Ca22chr{chrom}A_C_albicans_SC5314"

def process_fasta(input_file, output_file):
    """
    Process FASTA file and convert headers.
    
    Args:
        input_file (str): Path to input FASTA file
        output_file (str): Path to output FASTA file
    """
    try:
        with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
            for line in infile:
                if line.startswith('>'):
                    converted_header = convert_header(line)
                    outfile.write(converted_header + '\n')
                else:
                    outfile.write(line)
        print(f"Successfully converted headers. Output written to {output_file}")
        
    except FileNotFoundError:
        print(f"Error: Input file '{input_file}' not found.")
    except Exception as e:
        print(f"Error processing file: {e}")

def main():
    if len(sys.argv) != 3:
        print("Usage: python convert_fasta_headers.py <input_fasta> <output_fasta>")
        print("Example: python convert_fasta_headers.py input.fasta output.fasta")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    
    process_fasta(input_file, output_file)

if __name__ == "__main__":
    main()
