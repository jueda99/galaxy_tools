# This script generates the input files for MINTmap, given a tRNA FASTA file, a CSV file with anticodon loop positions, a genome FASTA file, and a tRNA BED file.
# It performs the following steps:
# 1. Generate all possible tRF sequences from the tRNA sequences, considering all possible
# 2. Identifying whether each tRF maps exclusively to the tRNA space or also to other regions of the genome, by mapping them to the genome and applying an exonic mask based on the tRNA annotations.
# 3. Write the lookup table and the tRF types file for MINTmap.
# 
# Bowtie is used for mapping, and the script assumes that bowtie is installed and available in the system PATH.
#
# run with:
#   python scripts/mintmap_input_generator.py "data/dmel-all-tRNA-r6.61.fasta" "data/anticodon_loop_positions.csv" "data/dmel-all-chromosome-r6.61.fasta" "data/tRNA_chr.bed" "output_test_branch/lookup_table.txt" "output_test_branch/tRF_types.txt"


import subprocess
from Bio import SeqIO
import pandas as pd
from pathlib import Path
import argparse # for readindg command line arguments
import os # for creating folders
import hashlib # for generating md5sum
from collections import defaultdict # for creating dict for serching tRFs and their excllusivity
import shutil # for removing folders

######## INPUT ########


parser = argparse.ArgumentParser(description="Generate input files for MINTmap")

parser.add_argument("trna_fasta", help="tRNA FASTA file")
parser.add_argument("anticodon_file", help="CSV file with anticodon loop positions")
parser.add_argument("genome_fasta", help="genome FASTA file")
parser.add_argument("trna_bed", help="tRNA BED file")

parser.add_argument("lookup_table", help="lookup table file")
parser.add_argument("tRF_types", help="tRF types file")

parser.add_argument("--threads", type=int, help="number of threads to use for bowtie mapping, default is 8", default=8)

parser.add_argument("--cleanup", action="store_true", help="Include this flag to clean up intermediate files")


args = parser.parse_args()


TRNA_FASTA = args.trna_fasta
ANTICODON_FILE = args.anticodon_file
GENOME_FASTA = args.genome_fasta
TRNA_BED = args.trna_bed

THREADS = args.threads
CLEANUP = args.cleanup

######## OUTPUT ########

LOOKUP = args.lookup_table
TYPES_OUT = args.tRF_types


######## INTERMEDIATE FILES ########

INTERMEDIATE_FOLDER="intermediate_files/"
if not os.path.exists(INTERMEDIATE_FOLDER):
    os.mkdir(INTERMEDIATE_FOLDER)

TRNA_OUT = INTERMEDIATE_FOLDER + "tRNA_sequences.fasta"
MASK_FILE = INTERMEDIATE_FOLDER + "mask.txt"
TRF_FASTA = INTERMEDIATE_FOLDER + "trfs.fasta"

BOWTIE_IDX_DIR= INTERMEDIATE_FOLDER + "bowtie_index/"
BOWTIE_IDX= BOWTIE_IDX_DIR + "genome_idx"


######## CLASSIFication of tRFs########
#
# Clasify tRFs based on the start/end positions of the tRFs on the tRNA sequence and the position of the anticodon loop (if known).
# Return the type of the tRF
#
def classify_trf(start, end, offset, ac_start, full_trf_len):

    # if we know the anticodon loop position
    if ac_start != -1:
        if start in [0, 1] and (ac_start - 2) <= end <= (ac_start + 1): # if tRF starts at the 1st or 2nd pos of 5'end of tRNA and ends in ac loop position -2 ~ +1
            if offset == 1:
                return "5'-half"
            elif offset == 0 and start == 0:
                return "5'-half"
            else:
                return "i-tRF"

        if (ac_start - 1) <= start <= (ac_start + 2) and end in [full_trf_len, full_trf_len-1, full_trf_len - 2, full_trf_len - 3]: # if tRF starts at ac loop position -1 ~ +2 and ends at the last or penultimate pos of the tRNA
            return "3'-half"

    # if we don't know the anticodon loop position and the tRF starts at the 1st or 2nd pos of the tRNA
    if start in [0, 1]:
        if offset == 1:     # if a nucleotide is added at the 5' end of the tRNA sequence, then the 2nd position is still the "biological 1st position"
            return "5'-tRF"
        elif offset == 0 and start == 0:    # if no nucleotide is added at the 5' end of the tRNA sequence, then only the 1st position can be considered as the "biological 1st position"
            return "5'-tRF"
        else:
            return "i-tRF"

    # if we don't know the anticodon loop position and the tRF ends at the last or penultimate pos of the tRNA
    if end in [full_trf_len, full_trf_len-1, full_trf_len - 2, full_trf_len - 3]: # we consider that a tRF can end at the last position of the tRNA (which is the "biological" 3' end) or at the penultimate position (which is the "biological" 3' end before CCA addition)
        return "3'-tRF"

    return "i-tRF"



######## GENERATE tRFs ########
#
# Generate all possible unique tRF sequences from the tRNA sequences, considering all possible post-transcriptional modifications (by adding A/T/C/G or nothing at the 5' end and CCA at the 3' end of the tRNA sequence) and classifying them with the classify_trf function.
# Fills the trf dict with tRF id as keys and their sequence with/without modifications and types as values.
# 
#
def generate_trfs(seq, ac_start, trf_dict, chr_strand_pos, seen_trf_seqs):

    variants = ["", "A", "T", "C", "G"]     # nucleotide added at the 5' end of the tRNA sequences to accommodate post-transcriptional modifications

    seq = seq + "CCA"   # add CCA, a tRNA maturation signal, at the 3' end of each tRNA sequence


    for prefix in variants:
        s = prefix + seq    # s = A/T/G/C/nothing + raw tRNA sequence + CCA
        offset = len(prefix)
        for i in range(len(s)-15):    # - 16 because the minimum length of tRFs is 16
            for l in range(16, 51):     # range starts at 16 because the minimum length of tRFs is 16

                start = i
                end = i + l

                if end > len(s): 
                    break

                s_e = str(start) + "_" + str(end)

                ac = ac_start + offset if ac_start != -1 else -1    # if anticodon position is known, ac_start contains the position, otherwise it contains "-1"

                t = classify_trf(start, end, offset, ac, len(s)) 


                full_trf = s[start:end]   # full tRF sequence, which can contain the post-transcriptional modifications (A/T/C/G at the 5' end and CCA at the 3' end)
                
                # get start and end position without prefix & CCA modifications
                clean_start = (start + 1) if (offset == 1 and start == 0) else start
                clean_end = min(end, len(s) - 3) 

                clean_fragment = s[clean_start:clean_end]

                # filtering 
                if full_trf in seen_trf_seqs: # if we have already seen this tRF sequence, we skip it to avoid duplicates in the output
                    continue
                
                seen_trf_seqs.add(full_trf)

                # Unique ID for each tRF
                prefix_id = "0" if offset == 0 else prefix
                s_e = f"{start}_{end}"
                tRFid = f"{prefix_id}_{chr_strand_pos}_{s_e}"
                
                trf_dict[tRFid] = [full_trf, clean_fragment, t]


    return



####### PARSE ANTICODON LOOP POSITIONS ########
def parse_anticodon_positions(ANTICODON_FILE):
    df = pd.read_csv(ANTICODON_FILE)
    anticodon_positions = {row["tRNA_Name"]: row["Loop_Start"] for _, row in df.iterrows()} # _ is the index of the row, we don't need it here, so we can use _ as a placeholder (equals to 'for index, row in df.iterrows()' but we only keep the row. For each row of the df take the row no matter what is its index)
    # this will do the same thing: anticodon_positions = dict(zip(df["tRNA_Name"], df["Loop_Start"]))
    return anticodon_positions



####### GET ANTICODON LOOP POSITIONS FOR MITOCHONDRIAL tRNAs #######
#
# Takes in charge of the mitochondrial tRNAs, then returns the start position of the anticodon loop if it can be found in the 16nt window around the middle of the tRNA sequence.
# Otherwise, returns -1 (which means the position of the anticodon loop is unknown)
#
def mt_tRNA_acloop_finder(trna_sequence, key, anticodon_dict):
    if key in anticodon_dict:
        mid = len(trna_sequence) // 2
        window = trna_sequence[max(0, mid-8):min(len(trna_sequence), mid+8)]  # generate a search window of 16nt in the center of the tRNA. use max() and min() just in case the tRNA sequence is shorter than 16 (should never happen)
        #print(f"trna_sequence\n{trna_sequence}\nkey = {key}\nwindow = {window}\nanticodon_dict[key] = {anticodon_dict[key]}\n")
        pos = window.find(key.split("-")[-1]) # search for the anticodon sequence and get position. key is in format: tRNA-AminoAcid-Codon
        #print(f"pos = {pos}\n")
        # if anticodon found in the window, calculate the start position in the whole tRNA sequence and return it; otherwise, return -1
        if pos != -1:
            ac_start_pos = max(0, mid-8) + pos
            return ac_start_pos
        else:
            return -1
    else:
        return -1


######## EXONIC MASK ########
#
# Build an exonic mask file for each chromosome sequence:
#   - initialize all positions to 0 (non-tRNA regions)
#   - set positions overlapping annotated tRNA exons to 1 using chromosome/strand/coordinates
#   - set positions corresponding to post-transcriptional additions (−1 and CCA) to 2
# Final mask encoding: 0 = background, 1 = tRNA exon, 2 = post-transcriptional addition
#
def create_mask():

    genome = list(SeqIO.parse(GENOME_FASTA, "fasta"))
    masks = {}
    # masks is a dict of the form 
    #   {chr1: [mask_forward_strand, mask_reverse_strand],
    #    chr2: [mask_forward_strand, mask_reverse_strand],
    #    ...}
    #

    # Initiate masks for each chromosome for both strands with 0
    for rec in genome:
        L = len(rec.seq)
        # for both 
        masks[rec.id] = [bytearray(b'0' * L), bytearray(b'0' * L)] # we use bytearray for memory efficiency, since the masks are very large
        '''
        masks[rec.id] = [["0"] * L]
        masks[rec.id].append(["0"] * L)
        '''

    with open(TRNA_BED) as f:
        next(f) # skip header
        for line in f:
            chrom, start, end, name, _, strand = line.strip().split()
            start, end = int(start)-1, int(end)-1 # convert to 0-based coordinates

            idx = [i for i, r in enumerate(genome) if r.id == chrom][0] # we suppose the genome fasta contains the chr information (which is stored in the "genome" SeqRecord object)

            if strand == "+": # modifying masks[chromosome][0] for the forward strand
                for i in range(start, end+1): # setting to 1 the positions overlapping annotated tRNA exons
                    masks[chrom][0][i] = 49     # ASCII for '1'
                for i in range(end+1, end+4): # setting to 2 the positions corresponding to post-transcriptional additions (CCA)
                    if i < len(masks[chrom][0]):
                        masks[chrom][0][i] = 50     # ASCII for '2'
                if start -1 >=0:
                    masks[chrom][0][start-1] = 50     # setting to 2 the position corresponding to one base upstream of tRNA

            else: # modifying masks[chromosome][1] for the reverse strand
                for i in range(start, end+1):
                    masks[chrom][1][i] = 49     # ASCII for '1'
                for i in range(start-3, start):
                    if i>= 0:
                        masks[chrom][1][i] = 50     # ASCII for '2'
                if end + 1 < len(masks[chrom][1]):
                    masks[chrom][1][end + 1] = 50     # ASCII for '2'
    return masks

    '''
    with open(MASK_FILE, "w") as f:
        for m in masks:
            two_line = f'>{m}_strand_+\n{"".join(masks[m][0])}\n>{m}_strand_-\n{"".join(masks[m][1])}\n'
            f.write(two_line)
    '''



######## BUILD BOWTIE INDEX ########
def build_index():
    Path(BOWTIE_IDX_DIR).mkdir(parents=True, exist_ok=True)
    subprocess.run(["bowtie-build", GENOME_FASTA, BOWTIE_IDX], check=True)



######## tRF MAPPING ########
#
# Bowtie mapping coommand to map the tRF sequences to the genome, allowing only perfect matches (no mismatch allowed, -v 0) and reporting all alignments (-a).
#
def map_trfs(fasta):
    return subprocess.run([
        "bowtie", "-v", "0", "-a",
        "--threads", str(THREADS),
        BOWTIE_IDX,
        "-f", fasta
    ], capture_output=True, text=True)



######## PARSE BOWTIE OUTPUT ########
# Parse the output of bowtie mapping.
# Return a dict :
#  {trf_id: [(chr, strand, position), (chr, strand, position), ...],
#   trf_id2: [(chr, strand, position), (chr, strand, position), ...],
#   ...}
#
def parse_hits(output):
    hits = {}
    for l in output.stdout.splitlines():
        # print(l)
        c = l.split("\t")
        trf_id = c[0]
        strand = c[1]
        chromos = c[2]
        pos = int(c[3]) # bowtie output is 0-based
        hits.setdefault(trf_id, []).append((chromos, strand, pos))
    return hits



######## APPLY EXONIC MASK ########
#
# Load the mask file into a dict of the form:
#       {chr1: [mask_forward_strand, mask_reverse_strand], chr2: [mask_forward_strand, mask_reverse_strand], ...}
#
def load_mask():
    mask_dict = {}
    with open(MASK_FILE) as f:
        lines = f.readlines() # convert in list so that we can use len() 
        for i in range(0,len(lines),2):
            chr_line = lines[i].strip(">").split("_strand_")[0]
            if chr_line not in mask_dict:
                mask_dict[chr_line] = [lines[i+1].strip()] 
            else:
                mask_dict[chr_line].append(lines[i+1].strip())


    return mask_dict



#
# Check if the tRF sequence is exclusive to the tRNA space by applying the mask to all its mapping positions on the genome.
# If there's one (or more) 0 in the tRF mapping segment, then that tF is not exclusive to the tRNA space, and return False. Otherwise, if all the mapping segments are composed of 1s and/or 2s, then the tRF is exclusive to the tRNA space, and return True.
#
def is_exclusive(seq, hits, mask):
    L = len(seq)

    for ref, strand, pos in hits:
        if strand == "+":
            segment = mask[ref][0][pos:pos+L]
        else:
            segment = mask[ref][1][pos:pos+L]
        if 48 in segment: # if there's a 0 in the segment, it means that the tRF doesn't map exclusively on the tRNA space
            #####
            #print(ref,strand, pos, "\n",seq, "\n", segment)
            #####
            return False

    return True



######## GENERATE MD5 SUM ########
#
# Generate the md5sum 
#
#
def md5(fname):
    hash_md5 = hashlib.md5()
    with open(fname, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash_md5.update(chunk)
    return hash_md5.hexdigest()




######## MAIN ########

def main():

    mask = create_mask() # tRNA mask: region of tRNA = 1, one base upstream and 3 bases downstream of tRNA = 2, other region = 0
    build_index()

    all_trfs = []
    types = {}
    
    trf_dict = {}    # dict of the form {tRF_id: tRFsequence, tRF sequence from unmodified tRNA, tRF_type}
    seen_trf_seqs = set()   # set to keep track of the tRF sequences to avoid duplicates

    # parse anticodon loop positions into a dict {tRNA_name: anticodon_loop_start_position}
    anticodon_dict = parse_anticodon_positions(ANTICODON_FILE) 

    # handle tRNA one by one to generate tRFs and classify them
    for rec in SeqIO.parse(TRNA_FASTA, "fasta"):
        seq = str(rec.seq)
        #print(f"rec = \n{rec.description}\n")
        
        # extract the tRNA name and for the name used in the anticodon position file
        #anticodon_key = "tRNA-"+"-".join(rec.description.split("name=")[1].split(";")[0].split(":")[-1].split("-")[:2])
        anticodon_key = "tRNA-"+"-".join(rec.description.split("_")[1].split("-")[:2])
        #print(anticodon_key)

        # extract tRNA starnd and postion for tRF_id definition: starnd_start_end
        chr_strand_pos = "_".join(rec.description.split('_')[-4:])
        #print(f"chr_strand_pos = {chr_strand_pos}")


        ### mitochondrial tRNAs
        #if rec.description.split("loc=")[1].startswith("mitochondrion_genome"):
        if rec.description.startswith("trnaMT"):
            ac_start_position = mt_tRNA_acloop_finder(seq, anticodon_key, anticodon_dict)
            trfs = generate_trfs(seq, ac_start_position, trf_dict, chr_strand_pos, seen_trf_seqs)

        ### nuclear tRNAs
        elif anticodon_key in anticodon_dict:
            #print(f"anticodon found for {anticodon_key} at position {anticodon_dict[anticodon_key]}")
            trfs = generate_trfs(seq, anticodon_dict[anticodon_key], trf_dict, chr_strand_pos, seen_trf_seqs)
        else:
            trfs = generate_trfs(seq, -1, trf_dict, chr_strand_pos, seen_trf_seqs) # if we don't know the anticodon loop position, we pass -1 to the generate_trfs function, which will then classify all tRFs as 5'-tRF, 3'-tRF or i-tRF based on their start and end positions only (without considering the anticodon loop position)

        #for s, t in trfs:   # s = tRF sequence, t = tRF type
        #    types[s] = t



    unique = list(types.keys()) # contains all unique tRF sequences



    unique_trfs = {}
    for trfid, values in trf_dict.items():
        clean_seq = values[1]  # get the "unmodified" sequence
    
        # If we haven't seen this biological sequence yet, keep it
        if clean_seq not in unique_trfs:
            unique_trfs[clean_seq] = {
                "id": [trfid],
                "raw_seq": values[0]
            }
        # otherwise, we only keep its id
        else:
            unique_trfs[clean_seq]["id"].append(trfid)

    with open(TRF_FASTA, "w") as f:    
        for clean_seq in unique_trfs:
            f.write(f">{('.').join(unique_trfs[clean_seq]['id'])}\n{clean_seq}\n")
        





    '''
    # write fasta of tRFs
    with open("old_fasta", "w") as f:
        for i, s in enumerate(unique):
            f.write(f">trf_{i}\n{s}\n")
    '''


    res = map_trfs(TRF_FASTA) # res contains the output of bowtie mapping of tRFs on the genome
    #print(res)
    hits = parse_hits(res)  # hits is a dict of mapping tRF IDs (e.g. "trf_0") to lists of (ref, pos) tuples where they map in the genome
    


    with open(TYPES_OUT, "w") as f:
        for i in trf_dict:
            f.write(f"{trf_dict[i][0]}\t{trf_dict[i][-1]}\n")
    

    with open(LOOKUP, "w") as out:

        # write md5sum on the header of the lookup file
        tRNA_md5sum = md5(TRNA_FASTA)
        trf_types_md5sum = md5(TYPES_OUT)
        out.write(f'#TRNASEQUENCES:{Path(TRNA_FASTA).name} MD5SUM:{tRNA_md5sum}\n')
        out.write(f'#OTHERANNOTATIONS:{Path(TYPES_OUT).name} MD5SUM:{trf_types_md5sum}\n')

        trf_to_hits_map = defaultdict(list)

        for h_key in hits:
            # Split the long compound key into individual tRF IDs
            # e.g., "0_X_..._0_16.A_X_..._0_17" -> ["0_X_..._0_16", "A_X_..._0_17"]
            individual_trf_ids = h_key.split('.')
            
            for t_key in individual_trf_ids:
                # Instantly check if this specific ID exists in your trf_dict (O(1) complexity)
                if t_key in trf_dict:
                    trf_to_hits_map[t_key].append(h_key)


        for i in trf_dict:
            # Retrieve all hit keys from the look-up map
            matching_hits_keys = trf_to_hits_map.get(i, [])
            if matching_hits_keys:
                # Evaluate exclusivity across ALL matching keys ()'is_exclusive' must be True for ALL matches to write 'Y')
                all_exclusive = True
                for h_key in matching_hits_keys:
                    if not is_exclusive(trf_dict[i][0], hits[h_key], mask):
                        all_exclusive = False
                        break  # As soon as one non-exclusive is found, stop checking 
                        
                out.write(f"{trf_dict[i][0]}\t{'Y' if all_exclusive else 'N'}\n")
            else:
                out.write(f"{trf_dict[i][0]}\tN\n")
    
    if CLEANUP:
        shutil.rmtree(INTERMEDIATE_FOLDER)  # remove the intermediate files folder and its contents

if __name__ == "__main__":
    main()
