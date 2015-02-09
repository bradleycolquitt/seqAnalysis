def format_fasta(fa):
    """Format fasta string for 50 characters per line"""
    max_ind = len(fa)/50
    fa1 = [fa[(x*50):(x*50)+50] + "\n" for x in range(0, max_ind)]
    fa1.append(fa[-1*(len(fa) - max_ind * 50) : -1])
    return fa1
