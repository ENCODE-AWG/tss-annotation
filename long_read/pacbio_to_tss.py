#!/usr/bin/python3
# Copyright 2020 California Institute of Technology
# Authors: Diane Trout
#
# BSD license with the additional requirement that
# Neither the name of the copyright holder nor the names of its
# contributors may be used to endorse or promote products derived from
# this software without specific prior written permission.
#
"""Generate TSS regions from an ENCODE PacBio transcript BAM file.
"""
from argparse import ArgumentParser

import pandas
import numpy
from typing import Union, List
from pysam import AlignmentFile, AlignedRead, idxstats
from collections import Counter, namedtuple
from heapq import heappush, heappop
import logging
from io import StringIO


DEFAULT_THRESHOLD = 20
DEFAULT_WINDOW = 20
MINUS_COLOR = "64,131,255"
PLUS_COLOR = "255,170,0"

logger = logging.getLogger("pacbio_to_tss")


tss_regions = namedtuple(
    "tss_regions",
    [
        "reference",
        "begin",
        "end",
        "is_reverse",
        "dense_begin",
        "dense_end",
        "summit",
        "expression",
    ],
)


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

    if args.chrom_info is not None:
        chrom_info = read_chrom_info(args.chrom_info, alignment.references)
    else:
        chrom_info = None

    reference_counts = get_counts_by_reference(args.bam_file)
    total_counts = reference_counts["mapped"].sum()
    logger.info("{} total reads".format(total_counts))

    threshold = args.expression_threshold
    if not args.raw_counts:
        raw_threshold = threshold * total_counts / 1_000_000
        threshold = int(numpy.ceil(raw_threshold))
        logger.info("Using {} as normalized threshold instead of {}".format(threshold, raw_threshold))

    regions, wigs = find_tss_peaks(
        alignment,
        threshold=threshold,
        window_size=args.window_size,
        use_tes=args.tes,
        chrom_info=chrom_info
    )

    if args.positive_bigwig is not None:
        update_wig(args.positive_bigwig, wigs[False], alignment)
    if args.negative_bigwig is not None:
        update_wig(args.negative_bigwig, wigs[True], alignment)
    wigs = None

    regions = pandas.DataFrame(regions)
    if not args.zero_bed:
        # We want [1,end] based beds. Internally everything is zero
        # based. So we need to add 1 to all the beginnings of
        # intervals to convert.
        regions["begin"] += 1
        regions["dense_begin"] += 1
        regions["summit"] += 1

    # Normalize expression to cpm
    if not args.raw_counts:
        regions["expression"] = regions["expression"] * 1_000_000 / total_counts

    # Calculate score field
    # score=int(math.log10(peakExp/maximum)*1000)
    peak_expression = regions["expression"].max()
    logger.info("Peak expression {}".format(peak_expression))
    regions["score"] = (
        numpy.ceil(regions["expression"] / peak_expression * 1000)
        .astype(int)
    )
    regions["strand"] = regions["is_reverse"].apply(lambda x: "-" if x else "+")
    regions["color"] = regions["is_reverse"].apply(
        lambda x: MINUS_COLOR if x else PLUS_COLOR
    )
    regions["name"] = ["{}{}".format(args.region_prefix, x) for x in regions.index]

    #          peak                           high density
    # chr, start, end, peak ID, score, strand, start, end, color, summit, expression
    columns = [
        "reference",
        "begin",
        "end",
        "name",
        "score",
        "strand",
        "dense_begin",
        "dense_end",
        "color",
        "summit",
        "expression",
    ]
    target = regions[columns]

    if args.output_file:
        target.to_csv(args.output_file, sep="\t", header=None, index=None)
    else:
        print(target)


def make_parser() -> ArgumentParser:
    parser = ArgumentParser(
        """Generate a bed file and optional wiggle files of likely tss regions

We assume that the input BAM file is zero based, and the values in the
bigwig will match the bam. By default the generated bed file will be
converted to 1 based coordinates.
"""
    )
    parser.add_argument(
        "--expression-threshold",
        type=float,
        default=DEFAULT_THRESHOLD,
        help="Require at least this many counts within a potential window",
    )
    parser.add_argument(
        "--window-size",
        type=int,
        default=DEFAULT_WINDOW,
        help="Expand window if the next read is within this many base pairs",
    )
    parser.add_argument("-i", "--bam-file", help="Transcript bam file", required=True)
    parser.add_argument("-o", "--output-file", help="TSS annotation bed file name")
    parser.add_argument("-n", "--negative-bigwig", help="minus strand bigwig file name")
    parser.add_argument("-p", "--positive-bigwig", help="plus strand bigwig file name")
    parser.add_argument(
        "--chrom-info",
        help="Use specified list of chromosomes. Helpful for creating UCSC compatible files"
    )
    parser.add_argument(
        "--tes",
        default=False,
        action="store_true",
        help="Use end site instead of start site",
    )
    parser.add_argument(
        "-0",
        "--zero-bed",
        default=False,
        action="store_true",
        help="Write out a zero based bed file.",
    )
    parser.add_argument(
        "-r",
        "--raw-counts",
        default=False,
        action="store_true",
        help="output raw counts instead of normalized to counts per million",
    )
    parser.add_argument(
        "--region-prefix", default="peak_", help="Set prefix for region names"
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


def get_counts_by_reference(filename: str) -> pandas.DataFrame:
    """Return table of mapped and unmapped read counts by reference name"""
    idx_buffer = StringIO(idxstats(str(filename)))
    names = ["reference_name", "length", "mapped", "unmapped"]
    counts = pandas.read_csv(idx_buffer, sep="\t", names=names)
    return counts


def find_tss_peaks(
    alignments: AlignmentFile,
    *,
    threshold: int = DEFAULT_THRESHOLD,
    window_size: int = DEFAULT_WINDOW,
    use_tes: bool = False,
    chrom_info: Union[None, List[str]] = None
) -> pandas.DataFrame:
    """Scan for tss peaks over all references in the BAM file.

    All input and output coordinates are assumed to be zero based.

    If wig_minus and wig_plus are pyBigWig files, the function will
    generate bigwigs of regions that pass the expression threshold.
    """
    bed_records = []
    wigs = {}

    if chrom_info is None:
        references = alignments.references
    else:
        references = chrom_info

    for name in sorted(references):
        records, chr_wigs = find_tss_peaks_on_reference(
            alignments, name, threshold=threshold, window_size=window_size, use_tes=use_tes
        )
        bed_records.extend(records)
        for strand in chr_wigs:
            wigs.setdefault(strand, {})[name] = chr_wigs[strand]

    return bed_records, wigs


def find_tss_peaks_on_reference(
    alignment: AlignmentFile,
    reference_name: str,
    *,
    threshold: int = DEFAULT_THRESHOLD,
    window_size: int = DEFAULT_WINDOW,
    use_tes: bool = False
) -> (List[tss_regions], {}):
    """Scan for tss peaks on a specific reference

    All input and output coordinates are assumed to be zero based.

    Returns a list of tss_regions and a dictionary of per base counts
    for tss_regions that pass the expression threshold for the plus
    and minus strand.
    (Basically what is needed to generate a bigwig file)
    """
    # Since we are using AlignedRead.is_reverse to determine what
    # strand a read is on, True means it is on the negative strand
    # while False means it's on the postive strand.
    window = {
        True: [],
        False: [],
    }
    wigs = {
        True: {},
        False: {},
    }
    bed_records = []

    for pos, read in fetch_next_read(alignment.fetch(reference_name), use_tes):
        strand = read.is_reverse

        # terminate a region
        if len(window[strand]) > 0 and pos > window[strand][-1] + window_size:
            summit, tss_count = find_most_frequent_tss(window[strand])
            if tss_count > threshold:
                bed_records.append(
                    calculate_tss_region(
                        reference_name, summit, tss_count, window[strand], strand
                    )
                )
                window_counts = Counter(sorted(window[strand]))
                for location in window_counts:
                    assert (
                        location not in wigs[strand]
                    ), "We should never see the same base twice {} {} {}".format(reference_name, strand, location)
                    wigs[strand][location] = window_counts[location]
            window[strand] = []
        # extend a region
        else:
            window[strand].append(pos)

    # We have to sort at the end because we built the positive and
    # negative strand windows separately so they can be a bit out of
    # order compared to what bedToBigBed wants
    bed_records = sorted(bed_records, key=lambda x: x.begin)
    logger.info("{} TSSes found on {}".format(len(bed_records), reference_name))
    logger.debug("{} minus strand wiggle on {}".format(len(wigs[True]), reference_name))
    logger.debug("{} plus strand wiggle on {}".format(len(wigs[False]), reference_name))
    return bed_records, wigs


def fetch_next_read(alignment, use_tes: bool) -> (int, AlignedRead):
    """Fetch reads in sorted order

    pysam returns reads ordered by reference_start, but we're
    interested in the transcription start sites, which for the reverse
    strand on the reference_end values.

    (or reverse if we're looking for the transcription end site (use_tes=True)
    """
    ordered_reads = []
    read_cache = {}

    for read in alignment:
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


def find_most_frequent_tss(region: List[int]) -> (int, int):
    """Find the summit and score the whole region

    Given a list of read ends, count them by base, and pick the maximum
    location, and return how many ends were in this region.
    """
    region_counts = Counter(region)
    sorted_counts = sorted(
        [(region_counts[location], location) for location in region_counts]
    )
    return sorted_counts[-1][1], len(region)


def calculate_tss_region(
    reference_name: str, summit: int, tss_count: int, region: List[int], is_reverse: bool
) -> tss_regions:
    median = numpy.median(region)
    sd = numpy.std(region)
    begin = region[0]
    end = region[-1] + 1
    assert begin < end, "We assume region is sorted in ascending order {} {}".format(begin, end)
    # Clamp the dense region to be inside the detected region. For
    # highly skewed regions the tail could fall outside the region
    dense_begin = max(int(numpy.floor(median - sd)), begin)
    dense_end = min(int(numpy.ceil(median + sd)), end)

    return tss_regions(
        reference_name,
        begin,
        end,
        is_reverse,
        dense_begin,
        dense_end,
        summit,
        tss_count,
    )


def add_wig_header_from_alignment(stream, alignment: AlignmentFile):
    """adds a BigWig header to the file using the BAM file's header"""
    stream.addHeader(list(sorted(zip(alignment.references, alignment.lengths))))


def update_wig(filename, wig, alignment):
    import pyBigWig

    stream = pyBigWig.open(filename, "w")
    add_wig_header_from_alignment(stream, alignment)

    for name in sorted(wig):
        logging.debug("Processing {} for wig".format(name))
        if len(wig[name]) > 0:
            starts = [int(x) for x in sorted(wig[name])]
            values = [float(wig[name][x]) for x in starts]
            stream.addEntries(name, starts, values=values, span=1)
    stream.close()


def read_chrom_info(filename, references):
    """Read chrom info file removing things not present in references
    """
    chrom_info = {}
    available = set(references)
    with open(filename, "rt") as instream:
        for line in instream:
            fields = line.strip().split("\t")
            if fields[0] in available:
                chrom_info[fields[0]] = int(fields[1])
            else:
                logger.warn(
                    "chrom_info name {} is not present in the bam file".format(fields[0])
                )

    return chrom_info


if __name__ == "__main__":
    main()
