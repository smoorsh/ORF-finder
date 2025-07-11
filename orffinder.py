#!/usr/bin/env python3

import re
from Bio import SeqIO
from Bio.Seq import Seq
import sys

##### Arguments ######
infile = sys.argv[1]
outfile = sys.argv[2]
min_length = sys.argv[3]

#### Proper input format:
# python3 orffinder.py /path/to/fasta/file.fa /path/to/output.gff minimum_bp_length_integer


###### Functions #######
def find_orfs(sequence, min_length):

    orfs = []
    sequence = sequence.upper()
    
    start_codon = 'ATG'
    stop_codons = set(['TGA', 'TAG', 'TAA'])  
    
    print(f"Analyzing sequence of length {len(sequence)}")
    
    # Check all 6 reading frames
    for frame in range(6):
        if frame < 3:
            seq = sequence[frame:]
            strand = '+'
            frame_num = frame + 1
            print(f"\nProcessing FORWARD frame {frame_num}")
        else:
            rev_seq = str(Seq(sequence).reverse_complement())
            seq = rev_seq[(frame-3):]
            strand = '-'
            frame_num = frame - 2
            print(f"\nProcessing REVERSE frame {frame_num}")
        
        # Process sequence codon by codon
        codon_positions = []
        for i in range(0, len(seq) - 2, 3):  # Move by 3 (codon by codon)
            codon = seq[i:i+3]
            codon_positions.append((i, codon))
        
        print(f"Total codons in frame: {len(codon_positions)}")
        
        # Find ORFs - each start gets exactly ONE ORF
        frame_orfs = 0
        i = 0
        
        while i < len(codon_positions) - 1:
            pos, codon = codon_positions[i]
            
            if codon == start_codon:
                start_pos = pos
                actual_start = start_pos + frame + 1 if strand == '+' else len(sequence) - (start_pos + frame)
                
                print(f"START CODON at codon_index {i}, seq_pos {pos}, actual_pos {actual_start}")
                
                # Look for the FIRST stop codon after this start
                found_stop = False
                j = i + 1  # Start from the NEXT codon
                
                while j < len(codon_positions):
                    stop_pos, stop_codon = codon_positions[j]
                    
                    if stop_codon in stop_codons:
                        # Found the FIRST stop codon
                        end_pos = stop_pos + 2  # Include the stop codon
                        actual_end = end_pos + frame + 1 if strand == '+' else len(sequence) - (end_pos + frame)
                        
                        orf_seq = seq[start_pos:stop_pos+3]  # Include stop codon
                        orf_length = len(orf_seq)
                        
                        print(f"FIRST STOP '{stop_codon}' at codon_index {j}, seq_pos {stop_pos}")
                        print(f"ORF length: {orf_length}bp")
                        
                        if orf_length >= int(min_length):
                            # Forward strand - CORRECT calculation
                            if strand == '+':
                                actual_start = i + frame + 1
                                actual_end = j + 2 + frame + 1
    
                            # Reverse strand - CORRECT calculation  
                            else:
                                actual_start = len(sequence) - (i + frame + 2)
                                actual_end = len(sequence) - (j + frame + 1)
                            
                            orfs.append((actual_start, actual_end, strand, frame_num, orf_seq))
                            frame_orfs += 1
                            print(f"ORF ADDED: {actual_start}-{actual_end} ({orf_length}bp)")
                        else:
                            print(f"ORF TOO SHORT: {orf_length} < {min_length}")
                        
                        found_stop = True
                        break  # CRITICAL: Stop searching after first stop
                    
                    j += 1
                
                if not found_stop:
                    print(f"No stop codon found for start at {actual_start}")
                
                print(f"Moving from codon_index {i} to {i+1}")
            
            i += 1  # Move to next codon position
        
        print(f"Frame {strand}{frame_num}: Added {frame_orfs} ORFs")
    
    print(f"\nFINAL BULLETPROOF RESULTS:")
    print(f"Total ORFs found: {len(orfs)}")
    
    # Ultra-strict duplicate detection
    start_check = {}
    duplicates_found = False
    
    for start, end, strand, frame, seq in orfs:
        key = (start, strand, frame)
        if key in start_check:
            duplicates_found = True
            print(f"CRITICAL BUG: Duplicate start {start} {strand} frame{frame}")
            print(f"   First ORF ends at: {start_check[key]}")
            print(f"   Duplicate ends at: {end}")
        else:
            start_check[key] = end
    
    if not duplicates_found:
        print(f"SUCCESS: No duplicate start positions!")
    
    return orfs

def write_gff_with_sequences(orfs, sequence_id, output_file):
    """
    Write ORFs to GFF file with embedded FASTA sequences
    """
    with open(output_file, 'w') as gff_file:
        # Write GFF header
        gff_file.write("##gff-version 3\n")
        gff_file.write(f"##sequence-region {sequence_id} 1 {len(sequence_id)}\n")
        
        # Write ORF features
        for i, (start, end, strand, frame, orf_seq) in enumerate(orfs, 1):
            # Calculate score based on length
            score = len(orf_seq)
            
            # Write GFF line
            gff_file.write(f"{sequence_id}\tORFFinder\tCDS\t{start}\t{end}\t{score}\t{strand}\t{frame-1}\t")
            gff_file.write(f"ID=ORF_{i:04d};Name=ORF_{i:04d};length={len(orf_seq)}\n")
        
        # Write FASTA section
        gff_file.write("##FASTA\n")
        
        # Write each ORF sequence
        for i, (start, end, strand, frame, orf_seq) in enumerate(orfs, 1):
            header = f">ORF_{i:04d} {sequence_id}:{start}-{end}({strand}) frame={frame} length={len(orf_seq)}"
            gff_file.write(f"{header}\n")
            
            # Write sequence in 80-character lines
            for j in range(0, len(orf_seq), 80):
                gff_file.write(f"{orf_seq[j:j+80]}\n")

def main():
    """Main program to find ORFs in FASTA file"""
    print(f"Starting ORF analysis...")
    print(f"Input file: {infile}")
    print(f"Output file: {outfile}")
    
    try:
        total_orfs = 0
        
        # Using BioPython
        for record in SeqIO.parse(infile, "fasta"):
            sequence_id = record.id
            sequence = str(record.seq)
            
            print(f"\nAnalyzing sequence: {sequence_id}")
            print(f"Sequence length: {len(sequence)} bp")
            
            # Find ORFs with algorithm
            orfs = find_orfs(sequence, min_length)
            
            print(f"\nFound {len(orfs)} ORFs >= {min_length}bp")
            
            if orfs:
                # Write to GFF file
                write_gff_with_sequences(orfs, sequence_id, outfile)
                total_orfs += len(orfs)
                
                # Print ORF summary - ONLY first 10 to avoid spam
                print(f"\nORF Summary for {sequence_id} (first 10):")
                for i, (start, end, strand, frame, orf_seq) in enumerate(orfs[:10], 1):
                    print(f"  ORF_{i:04d}: {start}-{end} {strand} frame{frame} ({len(orf_seq)}bp)")
                
                if len(orfs) > 10:
                    print(f"  ... and {len(orfs) - 10} more ORFs")
        
        print(f"\nAnalysis complete!")
        print(f"Total ORFs found: {total_orfs}")
        print(f"Results written to: {outfile}")
        
    except Exception as e:
        print(f"Error: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()