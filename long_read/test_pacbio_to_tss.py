#!/usr/bin/python3
from unittest import TestCase, main
from collections import namedtuple

from pacbio_to_tss import (
    find_most_frequent_tss,
    fetch_next_read,
    calculate_tss_region,
)

# This is a minimal mock of a pysam.AlignedRead
mock_read = namedtuple("MockRead", ["is_reverse", "query_name", "reference_start", "reference_end"])


class test_pacbio_to_tss(TestCase):
    def test_fetch_next_read_forward(self):
        data = [
            mock_read(False, "r1", 100, 1001),
            mock_read(False, "r2", 100, 990),
            mock_read(True, "r3", 101, 200),
            mock_read(False, "r4", 102, 1003),
        ]

        starts = []
        reads = []
        for start, read in fetch_next_read(data, use_tes=False):
            starts.append(start)
            reads.append(read)

        self.assertEqual(starts, [100, 100, 102, 200])

        # lets try again using the transcription end site
        starts = []
        reads = []
        for start, read in fetch_next_read(data, use_tes=True):
            starts.append(start)
            reads.append(read)

        self.assertEqual(starts, [101, 990, 1001, 1003])

    def test_fetch_next_read_reverse(self):
        data = [
            mock_read(True, "r1", 100, 1001),
            mock_read(True, "r2", 100, 990),
            mock_read(False, "r3",  101, 200),
            mock_read(True, "r4", 102, 1003),
        ]

        starts = []
        reads = []
        for start, read in fetch_next_read(data, use_tes=False):
            starts.append(start)
            reads.append(read)

        self.assertEqual(starts, [101, 990, 1001, 1003])

        # Lets try again using the end site
        starts = []
        reads = []
        for start, read in fetch_next_read(data, use_tes=True):
            starts.append(start)
            reads.append(read)

        self.assertEqual(starts, [100, 100, 102, 200])

    def test_find_most_frequent_tss_clear_cut(self):
        window = [100, 101, 101, 105, 106, 106, 106, 108]
        summit, tss_count = find_most_frequent_tss(window)
        self.assertEqual(summit, 106)
        self.assertEqual(tss_count, len(window))

    def test_find_most_frequent_tss_two_peaks(self):
        window = [100, 101, 101, 105, 106, 106, 108]
        summit, tss_count = find_most_frequent_tss(window)
        # It's sorted, so it will prefer the highest base?
        # should we do something different?
        self.assertEqual(summit, 106)
        self.assertEqual(tss_count, len(window))

    def test_calculate_tss_region_narrow(self):
        window = [100, 101, 101, 105, 106, 106, 106, 108]
        read = mock_read(True, "r1", 100, 110)

        is_reverse = False

        region = calculate_tss_region("chr1", 106, len(window), window, is_reverse, read)
        self.assertEqual(region.reference, "chr1")
        self.assertEqual(region.begin, 100)
        self.assertEqual(region.end, 109)
        self.assertEqual(region.is_reverse, is_reverse)
        self.assertEqual(region.dense_begin, 102)
        self.assertEqual(region.dense_end, 109)
        self.assertEqual(region.summit, 106)
        self.assertEqual(region.expression, len(window))
        self.assertEqual(region.query_name, read.query_name)

    def test_calculate_tss_region_wide(self):
        window = [100] * 10 + [140] + [150] * 10
        read = mock_read(True, "r10", 150, 160)
        is_reverse = False
        summit, tss_count = find_most_frequent_tss(window)

        region = calculate_tss_region("chr1", summit, tss_count, window, is_reverse, read)
        self.assertEqual(region.reference, "chr1")
        self.assertEqual(region.begin, 100)
        self.assertEqual(region.end, 151)
        self.assertEqual(region.is_reverse, is_reverse)
        self.assertEqual(region.dense_begin, 115)
        self.assertEqual(region.dense_end, 151)
        self.assertEqual(region.summit, 150)
        self.assertEqual(region.expression, len(window))
        self.assertEqual(region.query_name, read.query_name)


if __name__ == "__main__":
    main()
