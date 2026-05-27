# run with:
#   python ../mintmap_input_generator_v4.py "dmel-all-tRNA-r6.61.fasta" "anticodon_loop_positions.csv" "dmel-all-chromosome-r6.61.fasta" "tRNA_chr.bed" "lookup_table.txt' 'tRF_types.txt"
# ce script a été édité pour l'entrainemnt de galaxy. Il a étét modifié pour qu'il prenne en arguemnt les inputs (au lieu de les hardcoder)

import subprocess
from Bio import SeqIO
import pandas as pd
from pathlib import Path
import argparse # for readindg command line arguments
import shutil  # for removing intermediate folders
import hashlib # for generating md5sum
import os # for removing intermediate files

######## INPUT ########
parser = argparse.ArgumentParser(description="Generate input files for MINTmap")

parser.add_argument("trna_fasta", help="tRNA FASTA file")
parser.add_argument("anticodon_file", help="CSV file with anticodon loop positions")
parser.add_argument("genome_fasta", help="genome FASTA file")
parser.add_argument("trna_bed", help="tRNA BED file")
parser.add_argument("lookup_table", help="lookup table file")
parser.add_argument("tRF_types", help="tRF types file")
parser.add_argument("--threads", type=int, help="number of threads to use for bowtie mapping", default=8)


args = parser.parse_args()


TRNA_FASTA = args.trna_fasta
print(f"tRNA FASTA: {TRNA_FASTA}")
ANTICODON_FILE = args.anticodon_file
print(f"Anticodon loop positions file: {ANTICODON_FILE}")
GENOME_FASTA = args.genome_fasta
print(f"Genome FASTA: {GENOME_FASTA}")
TRNA_BED = args.trna_bed
print(f"tRNA BED file: {TRNA_BED}")

THREADS = args.threads

######## OUTPUT ########


LOOKUP = args.lookup_table
#TRNA_OUT = "tRNA_sequences.fasta"
TYPES_OUT = args.tRF_types

GENOME_TXT = "genome_space.txt"
MASK_FILE = "mask.txt"


######## OTHER PARAMETERS ########
BOWTIE_IDX_DIR= "bowtie_index/"
BOWTIE_IDX= BOWTIE_IDX_DIR + "genome_idx"


######## CLASSIFY ########

def classify_trf(start, end, ac_start, ac_end, length, trna_len):

    if ac_start != -1:
        if start in [0, 1] and (ac_start - 2) <= end <= (ac_start + 1):
            return "5'-half"

        if (ac_start - 1) <= start <= (ac_start + 2) and end in [trna_len, trna_len - 1]:
            return "3'-half"

    if start in [0, 1]:
        return "5'-tRF"

    if end in [trna_len, trna_len - 1]:
        return "3'-tRF"

    return "i-tRF"

######## GENERATE tRFs ########

def generate_trfs(seq, ac_start):

    variants = ["", "A", "T", "C", "G"]
    trfs = []

    seq = seq + "CCA"

    for prefix in variants:
        s = prefix + seq
        offset = len(prefix)

        for i in range(len(s) - 16):
            for l in range(16, 51):

                if i + l > len(s):
                    break

                start = i
                end = i + l

                ac = ac_start + offset if ac_start != -1 else -1

                t = classify_trf(start, end, ac, -1, l, len(s))

                trfs.append((s[start:end], t))

    return trfs

######## GENOME SPACE ########

def write_genome_space():
    with open(GENOME_TXT, "w") as out:
        for rec in SeqIO.parse(GENOME_FASTA, "fasta"):
            seq = str(rec.seq).upper()
            rev = str(rec.seq.reverse_complement()).upper()
            out.write(seq + "\n")
            out.write(rev + "\n")

######## EXONIC MASK ########

def create_mask():

    genome = list(SeqIO.parse(GENOME_FASTA, "fasta"))
    masks = {}

    for rec in genome:
        L = len(rec.seq)
        # for both strands
        masks[rec.id] = [["0"] * L]
        masks[rec.id].append(["0"] * L)

    with open(TRNA_BED) as f:
        next(f) # skip header
        for line in f:
            chrom, start, end, name, _, strand = line.strip().split()
            start, end = int(start), int(end)

            idx = [i for i, r in enumerate(genome) if r.id == chrom][0] # we suppose the genome fasta contains the chr information (which is stored in the "genome" SeqRecord object)

            if strand == "+":
                for i in range(start, end):
                    masks[chrom][0][i] = "1"
                for i in range(end, end+3):
                    if i < len(masks[chrom][0]):
                        masks[chrom][0][i] = "2"
                masks[chrom][0][start-1] = "2" # one base upstream of tRNA

            else:
                for i in range(start, end):
                    masks[chrom][1][i] = "1"
                for i in range(end, end+3):
                    if i < len(masks[chrom][1]):
                        masks[chrom][1][i] = "2"
                masks[chrom][1][start-1] = "2" # one base upstream of tRNA

    with open(MASK_FILE, "w") as f:
        for m in masks:
            two_line = f'>{m}_strand_+\n{"".join(masks[m][0])}\n>{m}_strand_-\n{"".join(masks[m][1])}\n'
            f.write(two_line)

######## BUILD BOWTIE INDEX ########

def build_index():
    Path(BOWTIE_IDX_DIR).mkdir(parents=True, exist_ok=True)
    subprocess.run(["bowtie-build", GENOME_FASTA, BOWTIE_IDX], check=True)

######## tRF MAPPING ########

def map_trfs(fasta):
    return subprocess.run([
        "bowtie", "-v", "0", "-a",
        "--threads", str(THREADS),
        BOWTIE_IDX,
        "-f", fasta
    ], capture_output=True, text=True)

######## PARSE BOWTIE OUTPUT ########

def parse_hits(output):
    hits = {}
    for l in output.stdout.splitlines():
        #print(l)
        c = l.split("\t")
        rid = c[0]
        strand = c[1]
        ref = c[2]
        pos = int(c[3])
        hits.setdefault(rid, []).append((ref, strand, pos))
    return hits

######## APPLY EXONIC MASK ########
#### changer en dict: [chr]: sequence mask
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
    #{m.split("\n")[0][1:]: m.split("\n")[1] for m in open(MASK_FILE).read().strip().split(">") if m}
    #[l.strip() for l in open(MASK_FILE)]

def is_exclusive(seq, hits, mask):

    L = len(seq)

    for ref, strand, pos in hits:
        if strand == "+":
            segment = mask[ref][0][pos:pos+L]
        else:
            segment = mask[ref][1][pos:pos+L]
        # print(ref,strand, pos, segment)
        if "0" in segment: # if there's a 0 in the segment, it means that the tRF doesn't map exclusively on the tRNA space
            return False

    return True

######## GENERATE MD5 SUM ########
def md5(fname):
    hash_md5 = hashlib.md5()
    with open(fname, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash_md5.update(chunk)
    return hash_md5.hexdigest()


######## MAIN ########

def main():

    write_genome_space()
    create_mask() # tRNA mask: region of tRNA = 1, one base upstream and 3 bases downstream of tRNA = 2, other region = 0
    build_index()

    all_trfs = []
    types = {}

    for rec in SeqIO.parse(TRNA_FASTA, "fasta"):
        seq = str(rec.seq)
        print(f"Processing tRNA: {rec.id}")
        trfs = generate_trfs(seq, -1)

        for s, t in trfs:   # s = tRF sequence, t = tRF type
            types[s] = t

    unique = list(types.keys()) # contains all unique tRF sequences

    # write fasta
    trfs_fasta = "trfs.fasta"
    with open(trfs_fasta, "w") as f:
        for i, s in enumerate(unique):
            f.write(f">trf_{i}\n{s}\n")

    res = map_trfs(trfs_fasta)
    #print(res)
    hits = parse_hits(res)  # hits is a dict mapping tRF IDs (e.g. "trf_0") to lists of (ref, pos) tuples where they map in the genome
    mask = load_mask()

    with open(TYPES_OUT, "w") as f:
        for s, t in types.items():
            f.write(f"{s}\t{t}\n")
    

    with open(LOOKUP, "w") as out:

        ####### add this in main script #######

        tRNA_md5sum = md5(TRNA_FASTA)
        trf_types_md5sum = md5(TYPES_OUT)
        out.write(f'#TRNASEQUENCES:{Path(TRNA_FASTA).name} MD5SUM:{tRNA_md5sum}\n')
        out.write(f'#OTHERANNOTATIONS:{Path(TYPES_OUT).name} MD5SUM:{trf_types_md5sum}\n')
        #########################################

        for i, s in enumerate(unique): # i = index of tRF sequence in unique, s = tRF sequence
            rid = f"trf_{i}"

            if rid not in hits: # if the tRF does not map on the genome; this is not supposed to happen, but just in case
                out.write(f"{s}\tN\n")
            else:
                ex = is_exclusive(s, hits[rid], mask)
                out.write(f"{s}\t{'Y' if ex else 'N'}\n")


    

    ####### add this in main script ? #######

    # remove intermediate files and folders
    shutil.rmtree(BOWTIE_IDX_DIR)
    os.remove(MASK_FILE)
    os.remove(GENOME_TXT)
    os.remove(trfs_fasta)

    #########################################

if __name__ == "__main__":
    main()