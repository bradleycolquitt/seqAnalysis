import pdb

def format_fasta(fa, start, width=50):
    """Format fasta string for width characters per line"""
    offset = width-start                   # To accomodate previously written lines
    max_ind = (len(fa) - offset)/width
    fa1 = [fa[:offset] + "\n"]
    fa1 = fa1 + [fa[(x*width) + offset : (x*width)+width + offset] + "\n" for x in range(0, (max_ind))]
    remainder = fa[max_ind * width + offset :]
    fa1.append(remainder)
    return (fa1, len(remainder))


def write_segment(seq, start, out, gap_size=1000):
    seq = seq.strip()
    ns = "N" * gap_size
    seq1 = seq + ns
    seq2, start = format_fasta(seq1, start)
    [out.write(x) for x in seq2]
    return start

complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
def reverse_complement(seq):
    bases = list(seq)
    bases = reversed([complement.get(base,base) for base in bases])
    bases = ''.join(bases)
    return bases
