#!/usr/bin/env python
import sys, random, itertools, argparse, os, logging
import HTSeq
from multiprocessing import Lock, JoinableQueue, Process


def process_chunk(q, l, fraction, out1, out2):
    while True:
        chunk = q.get()
        for read in chunk:
            if read and random.random() < fraction:
                read1 = HTSeq.SequenceWithQualities(*read[:3])
                read2 = HTSeq.SequenceWithQualities(*read[3:])
                l.acquire()
                read1.write_to_fastq_file(out1)
                read2.write_to_fastq_file(out2)
                out1.flush()
                out2.flush()
                l.release()
        q.task_done()

def grouper_nofill(n, iterable, iterable2=None):
    '''list(grouper_nofill(3, 'ABCDEFG')) --> [['A', 'B', 'C'], ['D', 'E', 'F'], ['G']]
    '''
    it=iter(iterable)
    if iterable2:
       it2=iter(iterable2)
    def take(paired):
        while 1:
            if paired:
                yield [(r1.seq,r1.name,r1.qualstr,r2.seq,r2.name,r2.qualstr) for r1,r2 in list(itertools.islice(itertools.izip(it, it2),n))]
            else:
                yield [(r.seq,r.name,r.qualstr) for r in list(itertools.islice(it,n))]
    return iter(take(iterable2 != None).next,[])

def run(fraction, input_fastq, output_dir, paired, procs, gzip_out, force):
    output_fastq = ["",""]
    fraction = float( fraction )

    in1 = iter(HTSeq.FastqReader( input_fastq[0] ))
    output_fastq[0] = "%s_%s_sampled.fastq" % (input_fastq[0].replace('.fastq','').replace('.gz',''),str(fraction))
    if output_dir:
        output_fastq[0] = os.path.join(os.path.abspath(output_dir),os.path.basename(output_fastq[0]))
    out1 = open( output_fastq[0], "w" )

    in2 = None
    out2 = None
    if paired:
        in2 = iter( HTSeq.FastqReader( input_fastq[1] ))
        output_fastq[1] = "%s_%s_sampled.fastq" % (input_fastq[1].replace('.fastq','').replace('.gz',''),str(fraction))
        if output_dir:
            output_fastq[1] = os.path.join(os.path.abspath(output_dir),os.path.basename(output_fastq[1]))
        out2 = open( output_fastq[1], "w" )

    lock = Lock()
    q = JoinableQueue()

    for i in range(procs):
        p = Process(target=process_chunk,args=(q, lock, fraction, out1, out2))
        p.daemon = True
        p.start()

    for chunk in grouper_nofill(1000, in1, in2):
        q.put(chunk)

    q.join()

    out1.close()
    if out2:
        out2.close()

    if gzip_out:
        import subprocess
        gzip_file = "%s.gz" % output_fastq[0]
        if os.path.exists(gzip_file) and not force:
            logging.warn("%s exists -- skipping" % gzip_file)
        else:
            subprocess.Popen(['pigz','-f','-p', str(procs), output_fastq[0]])
        if paired:
            gzip_file = "%s.gz" % output_fastq[1]
            if os.path.exists(gzip_file) and not force:
                logging.warn("%s exists -- skipping" % gzip_file)
            else:
                subprocess.Popen(['pigz','-f','-p', str(procs), output_fastq[1]])

if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog=__file__)

    parser.add_argument('-i', '--input', nargs='+', help="the input file[s]", required=True)
    parser.add_argument('-o', '--output-dir', help="the output folder. Default: same as input")
    parser.add_argument('-f', '--fraction', default="0.1", help="the fraction of read to subsample. Default: %(default)s")
    parser.add_argument('-p', '--processes', dest='procs', type=int, default=1, help="the number of processes. Default: %(default)s")
    parser.add_argument('-c', '--compress', action="store_true", default=False, help="gzip the output files. Default: %(default)s")
    parser.add_argument('--force', action="store_true", default=False, help="overwrite gzipped output files. Default: %(default)s")

    args = parser.parse_args()

    run(args.fraction, args.input, args.output_dir, len(args.input)==2, args.procs, args.compress, args.force)
