import argparse
import io
import unittest

from bin.get_targets_from_samfile import in_line

class TestProcessResults(unittest.TestCase):

    def test_in_line(self):
        target='19__RPOBgbs__RPOBgbs-4__19'
        line='@HD	VN:1.0	SO:unsorted\n'
        self.assertTrue(in_line(line, target))

        line='@SQ	SN:19__RPOBgbs__RPOBgbs-4__19	LN:33\n'
        self.assertTrue(in_line(line, target))

        line='@PG	ID:bowtie2	PN:bowtie2	VN:2.2.9	CL:"/opt/bowtie2-2.2.9/bowtie2-align-s --wrapper basic-0 -q --very-sensitive-local -a -x GBS_Res_Gene-DB_Final.fasta --passthrough -1 /tmp/29.inpipe1 -2 /tmp/29.inpipe2"'
        self.assertTrue(in_line(line, target))

        line='@SQ	SN:1__CAT__CAT-1__1	LN:100\n'
        self.assertFalse(in_line(line, target))

        line='HX6_25292:2:2113:29041:9150	145	19__RPOBgbs__RPOBgbs-4__19	9	255	25M126S	18__RPOBgbs__RPOBgbs-3__18	660	\
        TCCAAAAGCACCATATGTTGGTACTGGTATGGAGTATCAAGCAGCCCACGATTCAGGTGCAGCTGTGATTGCTAAACATGACGGTCGTGTTATTTTTTCAGATGCTGAAAAAGTTGAAGTGCGTCGCGAAGATGGTTCTCTTGATGTTTAT	\
        7JJFF<FAA<<JJJJJJJJJJJFFFFJJFJJJAJFAJJJJJJJJFFFFJJJJFJJJJJJJFFJJJJJJJJJJJJJJJFFJJFFJJJJJJAJJJJJJJFFFJFJJJJJJJJJJJJJJJFJF-JJJJJJJJJJJJJJJJJJJJJJJFJFF<AA	\
        AS:i:50	XN:i:0	XM:i:0	XO:i:0	XG:i:0	NM:i:0	MD:Z:25	YS:i:86	YT:Z:DP\n'
        self.assertTrue(in_line(line, target))

        line='HX6_25292:2:1216:15767:57389	137	2__ERMB__ERMB-1__2	1	255	80S71M	=	1	0	\
        GAAGGATTCTACAAGCGTACCTTGGATATTCACCGAACACTAGGGTTGCTCTTGCACACTCAAGTCTCGATTCAGCAATTGCTTAAGCTGCCAGCGGAATGCTTTCATCCTAAACCAAAAGTAAACAGTGTCTTAATAAAACTTACCCGCC	\
        AA<FFAJFJJJJJJJJJJJJJJJJJJJFJJJJJJJJJJJJJFJJJJJFJJFAFJJJJJJJJJJJJJJJJJJJFJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJFFJJJJJJJJJJJJJJJJJFFFJJFJAFFFAJFJJAJAJ<JJ-	\
        AS:i:142	XN:i:0	XM:i:0	XO:i:0	XG:i:0	NM:i:0	MD:Z:71	YT:Z:UP\n'
        self.assertFalse(in_line(line, target))
