#!/usr/bin/env python3
"""
AFChimera - AlphaFold for chimeric proteins with windowed MSA

AFChimera concatenates MSA files for chimeric protein structure prediction:
Standard windowed MSA: Tag-scaffold concatenation for N-terminus and C-terminus variants

Usage:
    python afchimera.py --scaffold-msa scaffold.a3m --tag-msa tag.a3m [options]
"""

import os
import argparse


def parse_a3m_file(filename):
    """Parse A3M file and return sequences dictionary"""
    sequences = {}
    with open(filename, 'r') as f:
        current_id = None
        for i, line in enumerate(f):
            if i == 0:  # Skip first line
                continue
            if line.startswith('>'):
                current_id = line.strip()
                sequences[current_id] = ''
            elif current_id is not None:
                sequences[current_id] += line.strip()
    return sequences


def get_terminus_tag(n_terminus):
    """Get terminus identifier string"""
    return "N" if n_terminus else "C"


def windowed_concatenation(file1, file2, output_file, n_terminus):
    """Windowed MSA concatenation"""
    sequences1 = parse_a3m_file(file1)  # scaffold
    sequences2 = parse_a3m_file(file2)  # tag

    with open(output_file, 'w') as out:
        all_ids = sorted(set(sequences1.keys()) | set(sequences2.keys()))
        
        # Write the key 101 first
        if '>101' in all_ids:
            out.write('>101\n')
            seq1 = sequences1.get('>101', '-' * len(next(iter(sequences1.values()))))
            seq2 = sequences2.get('>101', '-' * len(next(iter(sequences2.values()))))
            if n_terminus:
                out.write(f"{seq2}{seq1}\n")
            else:
                out.write(f"{seq1}{seq2}\n")
            all_ids.remove('>101')

        for seq_id in all_ids:
            out.write(f"{seq_id}\n")
            seq1 = sequences1.get(seq_id, '-' * len(next(iter(sequences1.values()))))
            seq2 = sequences2.get(seq_id, '-' * len(next(iter(sequences2.values()))))
            if n_terminus:
                out.write(f"{seq2}{seq1}\n")
            else:
                out.write(f"{seq1}{seq2}\n")


def run_concatenation(args):
    """Run MSA concatenation"""
    
    # Validate required inputs
    if not args.scaffold_msa:
        print("Error: --scaffold-msa is required")
        return
    if not args.tag_msa:
        print("Error: --tag-msa is required")
        return
    
    # Check if input files exist
    if not os.path.exists(args.scaffold_msa):
        print(f"Error: Scaffold MSA file not found: {args.scaffold_msa}")
        return
    if not os.path.exists(args.tag_msa):
        print(f"Error: Tag MSA file not found: {args.tag_msa}")
        return
    
    # Set output directory
    out_msas_folder = args.out_msas_folder or "./out_msas"
    os.makedirs(out_msas_folder, exist_ok=True)
    
    # Generate output file name from input files if not specified
    if args.output_file:
        output_base = args.output_file
    else:
        scaffold_name = os.path.splitext(os.path.basename(args.scaffold_msa))[0]
        tag_name = os.path.splitext(os.path.basename(args.tag_msa))[0]
        output_base = f"{scaffold_name}_{tag_name}"
    
    # Standard windowed concatenation mode
    print("Running windowed MSA concatenation...")
    
    processed_count = 0
    for n_terminus in [True, False]:
        terminus_tag = get_terminus_tag(n_terminus=n_terminus)
        print(f"Processing: n_terminus={n_terminus}")
        
        # Concatenate
        output_file = f'{out_msas_folder}/{terminus_tag}_{output_base}.a3m'
        windowed_concatenation(args.scaffold_msa, args.tag_msa, output_file, n_terminus)
        print(f"Created: {output_file}")
        processed_count += 1
    
    print(f"Processed {processed_count} windowed concatenations")


def create_parser():
    """Create command line argument parser"""
    parser = argparse.ArgumentParser(
        description="AFChimera - AlphaFold for chimeric proteins with windowed MSA",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Standard windowed concatenation (creates N-terminal and C-terminal variants)
  python afchimera.py --scaffold-msa scaffold.a3m --tag-msa tag.a3m
  
  # With custom output folder and filename
  python afchimera.py --scaffold-msa scaffold.a3m --tag-msa tag.a3m --out-msas-folder ./results --output-file my_chimera
        """
    )
    
    parser.add_argument('--scaffold-msa', required=True,
                       help='Path to the scaffold sequence MSA file (.a3m)')
    parser.add_argument('--tag-msa', required=True,
                       help='Path to the tag MSA file (.a3m)')
    parser.add_argument('--out-msas-folder', default='./out_msas',
                       help='Output MSAs folder (default: ./out_msas)')
    parser.add_argument('--output-file',
                       help='Base name for output files (auto-generated if not specified)')
    
    return parser


def main():
    """Main function"""
    parser = create_parser()
    args = parser.parse_args()
    
    try:
        run_concatenation(args)
        print("Successfully completed windowed concatenation!")
        
    except Exception as e:
        print(f"Error during execution: {e}")
        import traceback
        traceback.print_exc()
        return 1
    
    return 0


if __name__ == "__main__":
    exit(main()) 