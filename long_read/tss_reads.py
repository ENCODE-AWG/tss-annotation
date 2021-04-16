#!/usr/bin/python3
# Copyright 2020 California Institute of Technology
# Authors: Diane Trout
#
# BSD license with the additional requirement that
# Neither the name of the copyright holder nor the names of its
# contributors may be used to endorse or promote products derived from
# this software without specific prior written permission.
#
"""Given a TSS region file extract the reads that support the region
"""
from argparse import ArgumentParser

import pandas
import os
import sys
from pysam import AlignmentFile, AlignedRead
from heapq import heappush, heappop
import logging


logger = logging.getLogger("tss_reads")


def main(cmdline=None):
    parser = make_parser()
    args = parser.parse_args(cmdline)

    if args.debug:
        logging.basicConfig(level=logging.DEBUG)
    elif args.verbose:
        logging.basicConfig(level=logging.INFO)
    else:
        logging.basicConfig(level=logging.WARN)

    alignment = AlignmentFile(args.bam_file, "rb")
    if not alignment.has_index():
        logger.warning("Running without an index is not supported")

    #          peak                           high density
    # chr, start, end, peak ID, score, strand, start, end, color, summit, expression
    columns = [
        "reference",
        "begin",
        "end",
        "region_name",
        "score",
        "strand",
        "dense_begin",
        "dense_end",
        "color",
        "summit",
        "expression",
    ]
    column_dtypes = {
        "begin": int,
        "end": int,
        "dense_begin": int,
        "dense_end": int,
        "summit": int,
    }
    regions = pandas.read_csv(
        args.region_file,
        sep="\t",
        names=columns,
        dtype=column_dtypes,
        header=None,
        index_col=None,
    )

    if not args.zero_bed:
        # Bed is frequently stored as [1, end] and internally
        # bam is zero based, so we need to convert non-zero beds
        # to [0, end)
        regions["begin"] -= 1
        regions["dense_begin"] -= 1
        regions["summit"] -= 1

    if args.output_file is not None:
        outstream = open(args.output_file, "wt")
    else:
        outstream = sys.stdout

    try:
        for supporting_read in find_tss_supporting_reads(
            alignment, regions, args.tes, args.dense
        ):
            outstream.write("\t".join(supporting_read))
            outstream.write(os.linesep)
    finally:
        if args.output_file is not None:
            outstream.close()


def find_tss_supporting_reads(alignment, regions, use_tes=False, use_dense=False):
    """Find the reads whose TSS (or TES) contributed to the region

    alignment is a pysam AlignmentFile
    regions is a pandas DataFrame contining an ENCODE TSS group annotation file
    use_tes is a boolean flag to indicate if we should use the end site instead
    use_dense is a boolean flag to indicate if we should use the more narrow dense region

    :returns:
      iterator of peak name and read names
    """
    for line_number, row in regions.iterrows():
        if use_dense:
            begin = row.dense_begin
            end = row.dense_end
        else:
            begin = row.begin
            end = row.end

        for pos, read in fetch_next_read(
            alignment, reference=row.reference, begin=begin, end=end, use_tes=use_tes
        ):
            strand = '-' if read.is_reverse else '+'
            if pos >= begin and pos < end:
                yield (
                    row.reference,
                    str(read.reference_start),
                    str(read.reference_end),
                    read.query_name,
                    row.region_name,
                    strand,
                )


def make_parser() -> ArgumentParser:
    parser = ArgumentParser(
        """Generate a bed file and optional wiggle files of likely tss regions

We assume that the input BAM file is zero based, and the values in the
bigwig will match the bam. By default the generated bed file will be
converted to 1 based coordinates.
"""
    )
    parser.add_argument("-i", "--bam-file", help="Transcript bam file", required=True)
    parser.add_argument("-r", "--region-file", help="TSS region file", required=True)
    parser.add_argument("-o", "--output-file", help="TSS annotation bed file name")
    parser.add_argument(
        "--chrom-info",
        help="Use specified list of chromosomes. Helpful for creating UCSC compatible files",
    )
    parser.add_argument(
        "--tes",
        default=False,
        action="store_true",
        help="Use end site instead of start site",
    )
    parser.add_argument(
        "--dense",
        default=False,
        action="store_true",
        help="Use the more narrow dense region",
    )
    parser.add_argument(
        "-0",
        "--zero-bed",
        default=False,
        action="store_true",
        help="The input bed file is in zero-based coordinates instead of the more common one based.",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        default=False,
        action="store_true",
        help="print informative progress messages",
    )
    parser.add_argument(
        "--debug",
        default=False,
        action="store_true",
        help="report detailed log messages",
    )
    return parser


def fetch_next_read(
    alignment,
    *,
    reference: str = None,
    begin: int = None,
    end: int = None,
    use_tes: bool
) -> (int, AlignedRead):
    """Fetch reads in sorted order

    pysam returns reads ordered by reference_start, but we're
    interested in the transcription start sites, which for the reverse
    strand on the reference_end values.

    (or reverse if we're looking for the transcription end site (use_tes=True)
    """
    ordered_reads = []
    read_cache = {}

    for read in alignment.fetch(reference, begin, end):
        current = read.reference_start
        if read.is_reverse:
            pos = read.reference_end if not use_tes else read.reference_start
        else:
            pos = read.reference_start if not use_tes else read.reference_end

        # Sadly the ordered heap can't sort the AlignedSegment objects
        # So I have to store them elsewhere and only put entites into the
        # heap that have comparison operators defined.
        heappush(ordered_reads, (pos, read.query_name))
        read_cache[read.query_name] = read

        while len(ordered_reads) > 0 and current >= ordered_reads[0][0]:
            pos, query_name = heappop(ordered_reads)
            yield (pos, read_cache[query_name])
            del read_cache[query_name]

    while len(ordered_reads) > 0:
        pos, query_name = heappop(ordered_reads)
        yield (pos, read_cache[query_name])
        del read_cache[query_name]


if __name__ == "__main__":
    main()
