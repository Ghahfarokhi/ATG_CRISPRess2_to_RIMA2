#!/usr/bin/env python3

# atg_crispresso2rima.py
# This python script converts CRISPResso2 results into RIMA variant tables.
# Authors:  Mike Firth
#           Amir Taheri-Ghahfarokhi

import argparse 
import zipfile
import glob
import sys
import re
import os
import pandas as pd

pd.set_option('display.max_columns', None)

def revcomp(seq):
    """ Reverse complement of a sequence """
    return seq.translate(str.maketrans('ACGTacgtRYMKrymkVBHDvbhd', 'TGCAtgcaYRKMyrkmBVDHbvdh'))[::-1]

def collect_crispreeso_running_info(crispresso_folder):
    """Read the crispresso_info json file and extract running info"""
    alleles_frequency_files = sorted(glob.glob(crispresso_folder + "/**/Alleles_frequency_table.zip", recursive=True))
    
    if len(alleles_frequency_files) == 0 and crispresso_folder == "./":
        print_help()
        print("No CRISPResso folder found in the current working directory!")
        sys.exit(0)
    elif len(alleles_frequency_files) == 0:
        print(f"No CRISPResso folder found in {crispresso_folder}")
        sys.exit(1)

    df = pd.DataFrame({'alleles_frequency_files' : alleles_frequency_files})
    
    df['site_name'] = ''
    df['info_json_files'] = ''
    df['running_log_files'] = ''
    df['amplicon_seq'] = ''
    df['guide_seq'] = ''
    df['cut_pos'] = ''
    df['window'] = ''
    df['offset'] = ''

    df['info_json_files'] = df['alleles_frequency_files'].str.replace('Alleles_frequency_table.zip','CRISPResso2_info.json')
    df['running_log_files'] = df['alleles_frequency_files'].str.replace('Alleles_frequency_table.zip','CRISPResso_RUNNING_LOG.txt')

    for index, row in df.iterrows():

        json_file = row['info_json_files']
        running_log_file = row['running_log_files']

        if os.path.exists(json_file):

            amplicon_pattern = r'"amplicon_seq": "([^"]+)"'
            guide_pattern = r'"guide_seq": "([^"]+)"'

            with open(json_file, 'r') as file:
                data = file.read()

            amplicon = re.findall(amplicon_pattern, data)
            guide = re.findall(guide_pattern, data)
            df.loc[index, 'amplicon_seq'] = amplicon[0]
            df.loc[index, 'guide_seq'] = guide[0]

        elif os.path.exists(running_log_file):

            with open(running_log_file, 'r') as file:
                data = file.read()
            
            amplicon = re.search(r'--amplicon_seq (\S+)', data)
            guide = re.search(r'--guide_seq (\S+)', data)
            df.loc[index, 'amplicon_seq'] = amplicon.group(1)
            df.loc[index, 'guide_seq'] = guide.group(1)

        else:
            print(f"ERROR: Could't find either of CRISPResso2_info.json or CRISPResso_RUNNING_LOG.txt files for sample: {row['site_name']}")
            sys.exit(0)
    
    unmatch = []
    for index, row in df.iterrows():
        match = re.search("CRISPResso_on_([^/]+)", row['alleles_frequency_files'])
        if match:
            df.at[index, 'site_name'] = match.group(1)
        else:
            unmatch.append(row['alleles_frequency_files'])

    if len(unmatch) > 0:
        print("sample_name for the following CRISPResso folder(s) wasn't determined:\n")
        print("\n".join(unmatch))
    
    df.dropna(subset=['site_name'], inplace=True)
    
    if len(df.index) == 0:
        print("None of the folders have a regex pattern like 'CRISPResso_on_([^/]+)', from which the site name could be extracted.")
        print("Nothing else to be done!")
        sys.exit(1)

    for index, row in df.iterrows():
        amplicon_seq = row['amplicon_seq']
        guide_seq = row['guide_seq']
        amplicon_seq = amplicon_seq.upper()
        guide_seq = guide_seq.upper()
        
        if guide_seq in amplicon_seq:
            cut_pos = amplicon_seq.find(guide_seq) + len(guide_seq) - 3
        elif revcomp(guide_seq) in amplicon_seq:
            cut_pos = amplicon_seq.find(revcomp(guide_seq)) + 3
        else:
            print(f"Site: {row['site_name']} Error: Guide sequence ({guide_seq}) doesn't exist in the reference sequence ({amplicon_seq}).")
            continue

        df.at[index, 'cut_pos'] = cut_pos
        df.at[index, 'window'] = amplicon_seq[cut_pos-20:cut_pos+20]
        df.at[index, 'offset'] = cut_pos-20+1

    df.dropna(subset=['cut_pos'], inplace=True)
    if len(df.index) == 0:
        print("No site left to process after checking if the guide sequences exist in amplicon sequences!")
        print("Nothing else to be done!")
        sys.exit(1)
    
    # todo: Check if the site names are unique
    duplicate_sites_df = df[df.duplicated(['site_name'])]

    if not len(duplicate_sites_df.index) == 0:
        print("The following sites are found to have dulicative site_names and removed from further analysis!")
        print(duplicate_sites_df[['alleles_frequency_files', 'site_name']])
        df.drop_duplicates(subset=['site_name'])

    if len(df.index) == 0:
        print("No site left to process after droping duplicative site_names!")
        print("Nothing else to be done!")
        sys.exit(1)
    
    return df

def determine_MMEJ_deletion(ref, aligned, Read_Status, n_deleted, n_inserted, n_mutated, nreads, pc, guide):
    
    # This function determines if a deletion is MMEJ associated
    dashes = "-" * n_deleted

    # Find all locations on the aligned sequence that start with dashes
    positions = []
    index = aligned.find(dashes)
    
    while index != -1:
        positions.append(index)
        # Find the next occurrence of dashes, starting after the current one
        index = aligned.find(dashes, index + 1)
    # print("Positions: ", positions)
    
    # Choose dashes closest to the sgRNA
    # Determine if guide is on Fwd or Rev
    if ref.find(guide) == -1:
        sgRNA_position = ref.find(revcomp(guide))
    else:
        sgRNA_position = ref.find(guide)

    if sgRNA_position == -1 or len(positions) == 0:
        return 0

    # print("sgRNA pos: ", sgRNA_position)
    if len(positions) == 1:
        selected_dashes_position = positions[0]
    else:
        selected_dashes_position = min(positions, key=lambda x: abs(x - sgRNA_position))

    # Left-align the deletions
    for i in range(n_deleted):
        location = aligned.find(dashes, selected_dashes_position)
        if location - 1 > 0 and  location + n_deleted - 1 < len(ref):
            if aligned[ location - 1 ] == ref[ location + n_deleted - 1 ]:
                aligned = aligned[ :location - 1 ] + "-" + aligned[ location: ]
                aligned = aligned[ :location + n_deleted -1 ] + ref[ location + n_deleted - 1 ] + aligned[ location + n_deleted :]
                selected_dashes_position -= 1
    
    MMEJ_len = 0
    for i in range(n_deleted):
        location = aligned.find(dashes, selected_dashes_position)
        if location + i > 0 and location + n_deleted + i < len(ref):
            if ref[ location + i ] == aligned[ location + n_deleted + i]:
                MMEJ_len += 1

    return MMEJ_len

def process_alignment(site, options, ref, aligned, nreads, pc, combo):
    """ process one line of CRISPResso alignment into insertions, deletions and SNPs """

    # Loop through ref and aligned sequences assigning variants as we go
    s = ""
    n = len(ref)
    vars = []
    align_pos = 0
    ref_pos = 0
    while align_pos < n:
        if ref[align_pos] == aligned[align_pos]:
            c = " "
            align_pos += 1
            ref_pos += 1
        elif ref[align_pos] == "-":
            c = "+"
            insertion = ""
            while align_pos < n and ref[align_pos] == "-":
                insertion += aligned[align_pos]
                align_pos += 1
            vars.append({ "type":"Insertion", "from_pos":ref_pos, "to_pos":ref_pos, "ref":"-", "alt":insertion, "length":len(insertion) })
        elif aligned[align_pos] == "-":
            c = "-"
            deletion = ""
            deletion_pos = ref_pos
            while align_pos < n and aligned[align_pos] == "-":
                deletion += ref[align_pos]
                align_pos += 1
                ref_pos += 1
            vars.append({ "type":"Deletion", "from_pos":deletion_pos, "to_pos":ref_pos-1, "ref":deletion, "alt":"-", "length":ref_pos-1-deletion_pos+1 })
        else:
            c = "|"
            if not(options['indels_only']):
                vars.append({ "type":"SNV", "from_pos":ref_pos, "to_pos":ref_pos, "ref":ref[align_pos], "alt":aligned[align_pos], "length":1 })
            align_pos += 1
            ref_pos += 1
        s += c

    # Print alignment if debugging
    """
    print(ref)
    print(s)
    print(aligned)
    print("")
    """

    # Convert the variants to string keys for summing the number of reads
    keys = []
    min_from_pos = n
    max_to_pos = 1

    for var in vars:
        from_pos = site["offset"]+var["from_pos"]
        to_pos = site["offset"]+var["to_pos"]
        # Skip variants outside the specified window around the cut
        if options['window'] != 0 and (to_pos < site["cut_pos"]-options['window'] or from_pos > site["cut_pos"]+options['window']):
            continue
        key = ",".join([var["type"], var["ref"], var["alt"], str(from_pos), str(to_pos), str(var["length"]) ])
        keys.append(key)
        min_from_pos = min(min_from_pos, from_pos)
        max_to_pos = max(max_to_pos, to_pos)

    if len(keys) == 0:
        key = "WT"
        min_from_pos = 1
        max_to_pos = n
    else:
        key = "|".join(keys)

    # Add in the new combination
    try:
        combo[key]["nreads"] += nreads
        combo[key]["pc"] += pc
    except:
        combo[key] = { "nreads":nreads,  "pc":pc, "from":min_from_pos, "to":max_to_pos }

def process_site(df, site, options):
    site_name = site['name']
    alleles_frequency_file = site['alleles_frequency_file']
    guide = site['guide']
    print(f"\tProcessing allele frequency table: {alleles_frequency_file}")

    zip = zipfile.ZipFile(alleles_frequency_file)
    f = zip.namelist()[0]
    with zip.open(f) as opened_alleles_frequency_file:
        
        # If we're analysing the full Allele_frequencies ignore the offset
        if "around_sgRNA" not in alleles_frequency_file:
            site['offset'] = 1

        # Loop through the rows in the Allele_frequecies table
        line_no = 0
        total_reads = 0
        MMEJ_reads = 0
        del_01_bp = 0
        del_02_05_bp = 0
        del_06_10_bp = 0
        del_above_10_bp = 0
        ins_01_bp = 0
        ins_02_05_bp = 0
        ins_06_10_bp = 0
        ins_above_10_bp = 0
        indel_combos = 0
        header = None
        combo = {}
        alleles_frequency_file_with_MMEJ_annotaion_file = alleles_frequency_file.replace(".zip", "_MMEJ_annotated.txt")
        
        if os.path.exists(alleles_frequency_file_with_MMEJ_annotaion_file):
            os.remove(alleles_frequency_file_with_MMEJ_annotaion_file)

        with open(alleles_frequency_file_with_MMEJ_annotaion_file, "a") as MMEJ_file:

            for line in opened_alleles_frequency_file:
                line = line.decode().rstrip("\n")
                line_no += 1
                a = line.split("\t")
                if line_no == 1:
                    header = a
                    MMEJ_file.write("\t".join(header) + "\tMMEJ" + "\n")
                else:
                    row = {}
                    for i in range(len(header)):
                        row[ header[i] ] = a[i]
                    ref_name = row.get("Reference_Name", "Reference")
                    ref = row["Reference_Sequence"]
                    aligned = row["Aligned_Sequence"]
                    Read_Status = row["Read_Status"]
                    n_deleted = int(row["n_deleted"])
                    n_inserted = int(row["n_inserted"])
                    n_mutated = int(row["n_mutated"])
                    nreads = int(row["#Reads"])
                    pc = float(row["%Reads"])
                    row["MMEJ_deletion"] = 0
                    
                    if n_inserted == 0:
                        if n_deleted == 1:
                            del_01_bp += nreads
                        elif n_deleted < 6:
                            del_02_05_bp += nreads
                        elif n_deleted < 11:
                            del_06_10_bp += nreads
                        elif n_deleted > 10:
                            del_above_10_bp += nreads
                    elif n_deleted == 0: 
                        if n_inserted == 1:
                            ins_01_bp += nreads
                        elif n_inserted < 6:
                            ins_02_05_bp += nreads
                        elif n_inserted < 11:
                            ins_06_10_bp += nreads
                        elif n_inserted > 10:
                            ins_above_10_bp += nreads
                    elif n_deleted > 0 and n_inserted > 0: 
                        indel_combos += nreads

                    if n_deleted >= 2 :
                        row["MMEJ_deletion"] = determine_MMEJ_deletion (ref, aligned, Read_Status, n_deleted, n_inserted, n_mutated, nreads, pc, guide)
                    
                    if row["MMEJ_deletion"] >= 2:
                        MMEJ_reads += nreads

                    MMEJ_file.write("\t".join(map(str, list(row.values()))) + "\n")

                    if ref_name != "Reference": continue
                    process_alignment(site, options, ref, aligned, nreads, pc, combo)

                    total_reads += nreads
            
        RIMA_variants_df = pd.DataFrame(columns=['VariantNo','Position','Type','Length','Ref','Alt','Count'])
        total_variant_rows = 0
        mod_no = 0
        for key in sorted(combo.keys(), key=lambda s:-combo[s]["nreads"]):
            var = combo[key]
            if var["nreads"] >= options['min_reads'] and var["pc"] >= 100.0 * options['min_allele_freguency']:
                mod_no += 1

                if key == "WT":
                    continue
                for type  in key.split("|"):
                    pos, length, ref, alt = 0, 0, "", ""
                    type, ref, alt, start, end, length = type.split(",")
                    if type == "SNV": # RIMA doesn't need SNVs - so this reduces the number of rows
                        continue
                    new_row = [mod_no, start, type, length, ref, alt, var["nreads"]]
                    RIMA_variants_df.loc[len(RIMA_variants_df)] = new_row
                    
                    total_variant_rows += 1
    
    tsv_file_name = f"RIMA_{site['name']}-variants.tsv"
    tsv_file_path = f"{options['out']}/{tsv_file_name}"
    
    RIMA_variants_df.to_csv(tsv_file_path, sep="\t", index=False)

    row_index = df.index[df['site_name'] == site_name].tolist()[0]

    df.loc[row_index, 'mapped_reads'] = total_reads
    df.loc[row_index, 'mmej_reads'] = MMEJ_reads

    df.loc[row_index, 'del_01_bp'] = del_01_bp
    df.loc[row_index, 'del_02_05_bp'] = del_02_05_bp
    df.loc[row_index, 'del_06_10_bp'] = del_06_10_bp
    df.loc[row_index, 'del_above_10_bp'] = del_above_10_bp
    df.loc[row_index, 'ins_01_bp'] = ins_01_bp
    df.loc[row_index, 'ins_02_05_bp'] = ins_02_05_bp
    df.loc[row_index, 'ins_06_10_bp'] = ins_06_10_bp
    df.loc[row_index, 'ins_above_10_bp'] = ins_above_10_bp
    df.loc[row_index, 'indel_combos'] = indel_combos

    df.loc[row_index, 'total_variants'] = total_variant_rows
    df.loc[row_index, 'wt_count'] = combo["WT"]["nreads"] if "WT" in combo else 0
    df.loc[row_index, 'file_name'] = tsv_file_name
    df.loc[row_index, 'file_address'] = f"C:\RIMA\Raw\{tsv_file_name}"
    
    print(f"\tVariant file for {site['name']} saved as {tsv_file_path}")
    return df

def print_help():
    description = """
    This script recursively collects all CRISPResso's "Alleles_frequence_table.zip" files in path 
    and converts them into variant tables that can be analysed using RIMA.
    """
    usage = """\nUsage:\n
    For printing the help:\n
    python atg_crispresso2rima.py  [-h]\n
    For converting CRISPResso results into RIMA variant tables:\n
    python atg_crispresso2rima.py [--crispresso_folder PATH_TO_CRISPResso_DIR]\n
    Options
    --crispresso_folder     [string]        default="./"
    --indels_only           [True/Fasle]    default=False
    --min_reads             [integer]       default=0
    --min_allele_freguency  [float]         default=0.0
    --window                [integer]       default=20
    --out                   [string]        default="./RIMA"
    Example:
    python atg_crispresso2rima.py --crispresso_folder ./ --window 8 --out CRISPResso_2_RIMA
    """
    references = """\nReferences:
    RIMA v1: Taheri-Ghahfarokhi A, Taylor BJM, Nitsch R, Lundin A, Cavallo AL, Madeyski-Bengtson K, Karlsson F, Clausen M, Hicks R, Mayr LM, Bohlooly-Y M, Maresca M. Decoding non-random mutational signatures at Cas9 targeted sites. Nucleic Acids Res. 2018 Sep 19;46(16):8417-8434. doi: 10.1093/nar/gky653. PMID: 30032200; PMCID: PMC6144780.
    RIMA v2: Wimberger S, Akrap N, Firth M, Brengdahl J, Engberg S, Schwinn MK, Slater MR, Lundin A, Hsieh PP, Li S, Cerboni S, Sumner J, Bestas B, Schiffthaler B, Magnusson B, Di Castro S, Iyer P, Bohlooly-Y M, Machleidt T, Rees S, Engkvist O, Norris T, Cadogan EB, Forment JV, Šviković S, Akcakaya P, Taheri-Ghahfarokhi A, Maresca M. Simultaneous inhibition of DNA-PK and Polϴ improves integration efficiency and precision of genome editing. Nat Commun. 2023 Aug 14;14(1):4761. doi: 10.1038/s41467-023-40344-4. PMID: 37580318; PMCID: PMC10425386.
    \nCodes prepared by Mike Firth and Amir Taheri-Ghahfarokhi
    \nPlease report bugs to: Amir.Taheri.Ghahfarokhi@Gmail.com
    """
    print(description, usage, references)
    sys.exit(0)
    
def parse_arguments():
    if '-h' in sys.argv or '--help' in sys.argv:
        print_help()
        sys.exit(0)
    
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument("--crispresso_folder",      default="./",      type=str,     help="The parent folder of CRISPResso results.")
    parser.add_argument("--indels_only",            default=False,     type=bool,    help="Exclude SNPs (default %default)")
    parser.add_argument("--min_reads",              default=0,         type=int,     help="Min. reads to report (default %default)")
    parser.add_argument("--min_allele_freguency",   default=0.0,       type=float,   help="Min. allele frequency to report (default %default)")
    parser.add_argument("--window",                 default=20,        type=int,     help="Window size (default %default)")
    parser.add_argument("--out",                    default="./RIMA",  type=str,     help="Output directory (default %default)")

    args = parser.parse_args()

    options = {
        "crispresso_folder" : args.crispresso_folder,
        "indels_only" : args.indels_only,
        "min_reads" : args.min_reads,
        "min_allele_freguency" : args.min_allele_freguency,
        "window" : args.window,
        "out" : args.out
    }

    # Make the output folder if it already doesn't exists
    os.makedirs(args.out, exist_ok=True)

    return options, args

def main():

    options, args = parse_arguments()

    df = collect_crispreeso_running_info(args.crispresso_folder)
    
    sites = {}

    for index, row in df.iterrows():
        sites[df.loc[index, 'site_name']] = {
            "name": df.loc[index, 'site_name'],
            "alleles_frequency_file" : df.loc[index, 'alleles_frequency_files'],
            "reference": df.loc[index, 'amplicon_seq'],
            "guide": df.loc[index, 'guide_seq'],
            "cut_pos": df.loc[index, 'cut_pos'],
            "window": df.loc[index, 'window'],
            "offset": df.loc[index, 'offset']
        }
    
    df['mapped_reads'] = ''
    df['total_variants'] = ''
    df['wt_count'] = ''
    df['file_name'] = ''
    df['file_address'] = ''

    for data, site in sites.items():
        print(site['name'])
        df = process_site(df, site, options)
        
    experiment_sheet_df = df[[  'file_address', 
                                'file_name', 
                                'amplicon_seq', 
                                'guide_seq', 
                                'mapped_reads', 
                                'wt_count', 
                                'total_variants', 
                                'cut_pos', 
                                'mmej_reads',
                                'del_01_bp',
                                'del_02_05_bp',
                                'del_06_10_bp',
                                'del_above_10_bp',
                                'ins_01_bp',
                                'ins_02_05_bp',
                                'ins_06_10_bp',
                                'ins_above_10_bp',
                                'indel_combos'  ]]
    
    experiment_sheet_tsv_file = f"{options['out']}/experiment_sheet.tsv"
    experiment_sheet_df.to_csv(experiment_sheet_tsv_file, sep="\t", index=False)
    print(f"\nTab-delimited experiment_sheet file saved as {experiment_sheet_tsv_file}\nDone!\n")

if __name__ == '__main__':
    main()