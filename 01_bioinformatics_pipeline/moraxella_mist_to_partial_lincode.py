#! /usr/bin/env python3
import json
import sys
import os

# Moraxella mismatch cutoffs
MISMATCH_CUTOFFS = [1290, 832, 281, 85, 20, 8, 4, 2, 1, 0]

def determine_kept_positions(nb_diff) -> int:
    kept = 0
    for i, threshold in enumerate(MISMATCH_CUTOFFS):
        if nb_diff <= threshold:
            kept = i + 1
    return kept

def create_lin_code(lin_code_in: str, kept: int) -> str:
    parts = lin_code_in.split('_')
    total_positions = len(parts)
    kept = min(kept, total_positions)
    return '_'.join(parts[:kept] + (['*'] * (total_positions - kept)))

if __name__ == '__main__':
    if len(sys.argv) < 3:
        print('Usage: python3 moraxella_mist_to_partial_lincode.py <json_file1> [json_file2 ...] <profiles.tsv>', file=sys.stderr)
        sys.exit(1)

    # The last argument is the TSV database; everything before it are JSON files
    tsv_file = sys.argv[-1]
    json_files = sys.argv[1:-1]

    # 1. Load the TSV database into memory (do this ONCE for speed)
    db_profiles = []
    with open(tsv_file) as f:
        header = f.readline().strip('\n').split('\t')
        loci_columns = header[1:-1]
        total_loci = len(loci_columns)

        for line in f:
            parts = line.strip('\n').split('\t')
            if len(parts) < len(header):
                continue
            
            db_profiles.append({
                'cgst': parts[0],
                'lincode': parts[-1],
                'alleles': parts[1:-1]
            })

    # Print the table header
    print("Sample\tClosest_cgST\tMismatches\tLoci_Matched\tMatch_Percent\tBase_LINcode\tPartial_LINcode\tPositions_Assigned")

    # 2. Process each JSON file
    for json_file in json_files:
        sample_name = os.path.basename(json_file)
        
        try:
            with open(json_file) as handle:
                data_in = json.load(handle)

            # Extract sample alleles and their tags
            sample_alleles = {}
            sample_tags = {}
            if 'alleles' in data_in:
                for locus, details in data_in['alleles'].items():
                    sample_alleles[locus] = details.get('allele_str', '')
                    tags = details.get('tags', [])
                    sample_tags[locus] = tags[0] if tags else 'MISSING'
            else:
                print(f"{sample_name}\tERROR: No allele data found", file=sys.stderr)
                continue

            # Compare against the loaded database
            min_diff = float('inf')
            best_matching_st = None
            best_matching_lincode = None
            best_matched_loci = 0

            for profile in db_profiles:
                mismatches = 0
                matched_loci = 0
                
                for i, locus in enumerate(loci_columns):
                    prof_al = profile['alleles'][i]
                    samp_al = sample_alleles.get(locus, '')
                    tag = sample_tags.get(locus, 'MISSING')
                    
                    if tag == 'INDEL' or tag == 'MISSING':
                        matched_loci += 1
                    elif tag == 'NOVEL':
                        mismatches += 1
                    elif tag == 'EXACT':
                        if prof_al == 'N':
                            matched_loci += 1
                        elif prof_al == samp_al:
                            matched_loci += 1
                        else:
                            mismatches += 1
                
                # Keep track of the closest profile
                if mismatches < min_diff:
                    min_diff = mismatches
                    best_matching_st = profile['cgst']
                    best_matching_lincode = profile['lincode']
                    best_matched_loci = matched_loci
                    
                    # Stop early if we find a perfect match
                    if mismatches == 0:
                        break

            if best_matching_st is None:
                print(f"{sample_name}\tERROR: Could not evaluate against DB", file=sys.stderr)
                continue

            # Calculate and format the output row
            match_percentage = (best_matched_loci / total_loci) * 100 if total_loci > 0 else 0
            kept_positions = determine_kept_positions(min_diff)
            partial_lin_code = create_lin_code(best_matching_lincode, kept_positions)
            
            # Print row data separated by tabs
            print(f"{sample_name}\t{best_matching_st}\t{min_diff}\t{best_matched_loci}/{total_loci}\t{match_percentage:.1f}\t{best_matching_lincode}\t{partial_lin_code}\t{kept_positions}")

        except Exception as e:
            # If one JSON is broken, log the error but don't stop the loop for the rest of the files
            print(f"{sample_name}\tFAILED: {str(e)}", file=sys.stderr)
