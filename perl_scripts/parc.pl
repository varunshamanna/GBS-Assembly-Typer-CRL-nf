#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper;
use Getopt::Std;
use File::Copy qw(copy);
use Env;

my %Res_Targets = (
    "ERM" => "neg",
    "LNUB" => "neg",
    "LSA" => "neg",
    "MEF" => "neg",
    "TET" => "neg",
    "CAT" => "neg",
    "GYRA" => "neg",
    "PARC" => "neg",
    "23S1" => "neg",
    "23S3" => "neg",
    "RPOB1" => "neg",
    "RPOB2" => "neg",
    "RPOB3" => "neg",
    "RPOB4" => "neg",
    );

my %drugRes_Col = (
    "TET" => "neg",
    "EC" => "neg",
    "FQ" => "neg",
    "OTHER" => "neg",
    );

my @Bin_Res_arr = (0) x 19;

sub extractFastaByID {
    my ($lookup, $reference) = @_;
    open my $fh, "<", $reference or die $!;
    local $/ = "\n>";  # read by FASTA record

    my $output;
    while (my $seq = <$fh>) {
        chomp $seq;
        #print "while seq:\n$seq\n";
        my ($id) = $seq =~ /^>*(\S+)/;  # parse ID as first word in FASTA header
        if ($id eq $lookup) {
            $seq =~ s/^>*.+\n//;  # remove FASTA header
            $seq =~ s/\n//g;  # remove endlines
            #print ">$id\n";
            #print "$seq\n";
            $output = ">$id\n$seq\n";
            #$output = $seq;
            last;
        }
    }
    return $output;
}

sub freebayes_prior_fix {
    my ($bamFile, $refFile, $target) = @_;
    (my $samFile = $bamFile) =~ s/\.bam/\.sam/g;
    system("samtools view -h $bamFile > $samFile");
    system("cat $samFile | grep -E \"^\@HD|^\@SQ.*$target|^\@PG\" > CHECK_target_seq.sam");
    system("awk -F'\t' '\$3 == \"$target\" {print \$0}' $samFile >> CHECK_target_seq.sam");
    system("samtools view -bS CHECK_target_seq.sam > CHECK_target_seq.bam");
    system("samtools index CHECK_target_seq.bam CHECK_target_seq.bai");
    my $REF_seq = extractFastaByID("$target","$refFile");
    open(my $rf,'>',"CHECK_target_ref.fna");
    print $rf "$REF_seq\n";
    close $rf;
    system("freebayes -q 20 -p 1 -f CHECK_target_ref.fna CHECK_target_seq.bam -v CHECK_target_seq.vcf");
    system("bgzip CHECK_target_seq.vcf");
    system("tabix -p vcf CHECK_target_seq.vcf.gz");
    my $extractSeq = `echo "$REF_seq" | vcf-consensus CHECK_target_seq.vcf.gz`;
    chomp($extractSeq);
    #print "$target-----------------------------------\n";
    #system("cat CHECK_target_seq.sam");
    #system("zcat CHECK_target_seq.vcf.gz");
    #print "reference seq:\n$REF_seq\n";
    #print "extracted Seq:\n$extractSeq\n";
    #print "$target-----------------------------------\n";
    system("rm CHECK_target*");
    return $extractSeq;
}

sub sixFrame_Translate {
    my ($seq_input,$opt_f) = @_;

    sub codon2aa{
	my($codon)=@_;
        $codon=uc $codon;
        my(%g)=('TCA'=>'S','TCC'=>'S','TCG'=>'S','TCT'=>'S','TTC'=>'F','TTT'=>'F','TTA'=>'L','TTG'=>'L','TAC'=>'Y','TAT'=>'Y','TAA'=>'*','TAG'=>'*','TGC'=>'C','TGT'=>'C','TGA'=>'*','TGG'=>'W','CTA'=>'L','CTC'=>'L','CTG'=>'L','CTT'=>'L','CCA'=>'P','CCC'=>'P','CCG'=>'P','CCT'=>'P','CAC'=>'H','CAT'=>'H','CAA'=>'Q','CAG'=>'Q','CGA'=>'R','CGC'=>'R','CGG'=>'R','CGT'=>'R','ATA'=>'I','ATC'=>'I','ATT'=>'I','ATG'=>'M','ACA'=>'T','ACC'=>'T','ACG'=>'T','ACT'=>'T','AAC'=>'N','AAT'=>'N','AAA'=>'K','AAG'=>'K','AGC'=>'S','AGT'=>'S','AGA'=>'R','AGG'=>'R','GTA'=>'V','GTC'=>'V','GTG'=>'V','GTT'=>'V','GCA'=>'A','GCC'=>'A','GCG'=>'A','GCT'=>'A','GAC'=>'D','GAT'=>'D','GAA'=>'E','GAG'=>'E','GGA'=>'G','GGC'=>'G','GGG'=>'G','GGT'=>'G');

        if(exists $g{$codon}){return $g{$codon};}
        elsif($codon=~/GC./i){return 'A';}
        elsif($codon=~/GG./i){return 'G';}
        elsif($codon=~/CC./i){return 'P';}
        elsif($codon=~/AC./i){return 'T';}
        elsif($codon=~/GT./i){return 'V';}
        elsif($codon=~/CG./i){return 'R';}
        elsif($codon=~/TC./i){return 'S';}
        else {
            return('x');
            print "Bad codon \"$codon\"!!\n";
        }
    }

    (my $DNAheader,my @DANseq)=split(/\n/,$seq_input);
    chomp $DNAheader;
    $DNAheader=~s/\s+$//g;
    my $DNAseq = join( '',@DANseq);
    $DNAseq =~ s/\s//g;
    $DNAheader=~s/>//g;
    $DNAseq=~s/>//g;
    my$DNA_length=length$DNAseq;
    #print "\nSeq:$DNAheader\t:$DNA_length nt\n\n";
    my $DNArevSeq = reverse($DNAseq);
    $DNArevSeq=~tr/ATGCatgc/TACGtacg/;
    #print "\nThe original DNA sequence is:\n$DNAseq \nThe reverse of DNA sequence is:\n$DNArevSeq\n";
    my @protein='';
    my @dna='';
    my $codon1;

    if ($opt_f == 1) {
        for(my $i=0;$i<(length($DNAseq)-2);$i+=3){
            $codon1=substr($DNAseq,$i,3);
            $protein[1].= codon2aa($codon1);
            #$dna[1].=codon2nt($codon1);
            }
        }
    if ($opt_f == 2) {
            my $codon2;
            for(my $i=1;$i<(length($DNAseq)-2);$i+=3){
                $codon2=substr($DNAseq,$i,3);
                $protein[2].= codon2aa($codon2);
                #$dna[2].=codon2nt($codon2);
                }
    }
    if ($opt_f == 3) {
            my $codon3;
            for(my $i=2;$i<(length($DNAseq)-2);$i+=3){
                $codon3=substr($DNAseq,$i,3);
                $protein[3].= codon2aa($codon3);
                #$dna[3].=codon2nt($codon3);
                }
    }
    if ($opt_f == 4) {
            my $codon4;
            for(my $i=0;$i<(length($DNArevSeq)-2);$i+=3){
                $codon4=substr($DNArevSeq,$i,3);
                $protein[4].= codon2aa($codon4);
                #$dna[4].=codon2nt($codon4);
                }
    }
    if ($opt_f == 5) {
            my $codon5;
            for(my $i=1;$i<(length($DNArevSeq)-2);$i+=3){
                $codon5=substr($DNArevSeq,$i,3);
                $protein[5].= codon2aa($codon5);
                #$dna[5].=codon2nt($codon5);
                }
    }
    if ($opt_f == 6) {
            my $codon6;
                for(my $i=2;$i<(length($DNArevSeq)-2);$i+=3){
                    $codon6=substr($DNArevSeq,$i,3);
                    $protein[6].= codon2aa($codon6);
                    #$dna[6].=codon2nt($codon6);
                    }
    }
#print "translate result\n$protein[$opt_f]\n";
return $protein[$opt_f];
}



$Res_Targets{"PARC"} = "pos";
my $RES_bam = "../test_data/RES_26189_8#5__26189_8#5.GBS_Res_Gene-DB_Final.sorted.bam";
my $res_DB = "../db/0.0.1/GBS_resTyper_Gene-DB/GBS_Res_Gene-DB_Final.fasta";

if ($Res_Targets{"PARC"} eq "pos") {
    my @PARC_output;
    my $PARC_seq;
    my $PARC_aaSeq;
    #my @miscR_value = split(':',$miscR_Type{"PARCGBS-1"});
    #my $PARC_seq = extractFastaByID("7__PARCGBS__PARCGBS-1__7","TEMP_miscR_consensus.fna");
    #my $PARC_seq = freebayes_prior_fix($RES_bam, $res_DB,"7__PARCGBS__PARCGBS-1__7");
    #my $PARC_aaSeq = sixFrame_Translate($PARC_seq,1);
    my $PARC_aaSeq = "HHHGDSSIYDAMVRMSS";
    my $PARC_aaRef = "HPHGDSSIYDAMVRMSQ";
    print "PARCgbs sequence: $PARC_seq || $PARC_aaSeq\n";
    if ($PARC_aaSeq ne "HPHGDSSIYDAMVRMSQ") {
        my $mask = $PARC_aaSeq ^ $PARC_aaRef;
        my @seq_diffs;
        while ($mask =~ /[^\0]/g) {
            print substr($PARC_aaRef,$-[0],1), ' ', substr($PARC_aaSeq,$-[0],1), ' ', $-[0], "\n";
            #my $diff_element = "pos".($-[0]+1).":".substr($RPOB_aaRef,$-[0],1)."->".substr($RPOB_aaSeq,$-[0],1);
            my $diff_element = substr($PARC_aaRef,$-[0],1).($-[0]+1).substr($PARC_aaSeq,$-[0],1);
            push(@seq_diffs,$diff_element);
        }
        print "PARC seq: $PARC_seq\n";
        my $diff_output = join(',',@seq_diffs);
        print "$diff_output\n";
        my $bin_out = join(':',@seq_diffs);
        print "$bin_out\n";
        $Bin_Res_arr[14] = $bin_out;
        my $PARC_out = "PARC-".$diff_output;
        print "$PARC_out\n";
        if ($drugRes_Col{"FQ"} eq "neg") {
            $drugRes_Col{"FQ"} = $PARC_out;
        } else {
            my $new_val = $drugRes_Col{"FQ"}.":".$PARC_out;
            $drugRes_Col{"FQ"} = $new_val;
        }
    }
}
