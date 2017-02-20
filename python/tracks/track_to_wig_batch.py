import sys, os
import argparse

import track_to_wig as tw
from multiprocessing import Pool

def worker(track_name, gdb=None, assembly=None):
    #print i
    #print in_files[i]
    #print out_files[i]
    #print assembly
    #tracks = [ in_files[i] + "-" + p for p in suffixes ]
    #gdb.delete_track()
    #print tracks
    tw.run(track_name, gdb, assembly)
    #combine_tracks.create_combined_tracks(out_files[i], tracks, assembly,
    #                       np.dtype(dtype))

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', dest='in_dir')
    parser.add_argument('-c', dest='ncores', type=int, help='number of threads to use', default=4)
    parser.add_argument("--dtype", metavar="", action="store",
                        choices=("uint8", "uint16"), default="uint8",
                        help="datatype of combined track")
    parser.add_argument('--assembly', help="assembly to use", default=None)
    args = parser.parse_args()

    gdb = genome.db.GenomeDB()
    all_files = gdb.list_tracks(subdir=args.in_dir, recursive=False)

    # trim_files = list(set([ re.sub(pattern, '', f) for f in all_files ]))
    # base_files = [os.path.basename(f) for f in trim_files]

    # out_files = ["/".join([args.in_dir, 'combined', f]) for f in base_files]
    # #pdb.set_trace()


    pool = Pool(processes=args.ncores)
    partial_worker = partial(worker,
                             gdb=gdb,
                             assembly=args.assembly)
    #partial_worker(0)
    pool.map(partial_worker, all_files, 1)
    #pool.close()
    # for i in xrange(len(out_files)):
    #    print out_files[i]
    #    tracks = [ trim_files[i] + "-" + p for p in suffixes]
    #    combine_tracks.create_combined_tracks(out_files[i], tracks, args.assembly,
    #                       np.dtype(args.dtype))


if __name__ == '__main__':
    main()
