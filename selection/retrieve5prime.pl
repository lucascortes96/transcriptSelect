#Perl Imports 
use Bio::EnsEMBL::Registry;
my $registry = 'Bio::EnsEMBL::Registry';
use Data::Dumper;
use Getopt::Long;

## Grab database we want to use 
$registry->load_registry_from_db(
    -host => 'ensembldb.ensembl.org',  # alternatively 'useastdb.ensembl.org'
    -user => 'anonymous'
);

##Grabbing genes and transcripts for protein coding genes 
my $gene_adaptor = $registry->get_adaptor('Human', 'Core', 'Gene');
my $transcript_adaptor = $registry->get_adaptor('Human', 'Core', 'Transcript');
my $genes = $gene_adaptor->fetch_all_by_biotype('protein_coding');


# Initialize two empty hashes to store the transcripts and geneIDs
my %mane_transcripts;
my %non_mane_transcripts;

##Looping through all genes in the human genome 
foreach my $gene (@{$genes}) {
    my $transcripts = $gene->get_all_Transcripts();
    #my $mane_transcript = $transcript_adaptor->fetch_MANE_by_gene_stable_id($gene->stable_id());
    #print $mane_transcript, "/n";
    foreach my $transcript (@{$transcripts}) {
        my $attributes = $transcript->get_all_Attributes('MANE_Select');
        if (@{$attributes}) {
            # Store the MANE Select transcript ID in the mane_transcripts hash
            push @{$mane_transcripts{$gene->stable_id}}, $transcript->stable_id;
        }
        else { 
            # Store the non-MANE Select transcript ID in the non_mane_transcripts hash
            push @{$non_mane_transcripts{$gene->stable_id}}, $transcript->stable_id;
        }
        #my $mane_transcript = $transcript_adaptor->fetch_by_stable_id($transcript->stable_id());
        #my $mane_gene_stable_id = $mane_transcript->get_Gene()->stable_id();
        #print $mane_transcript->stable_id(), ", ", $mane_gene_stable_id, ",", $transcript->stable_id(), "\n";
        #print $gene->stable_id(), "\n";
        #print "Transcript Stable ID: ", $transcript->stable_id(), "\n";
        #if ($mane_transcript && $mane_gene_stable_id eq $gene->stable_id()) {
            #print "match!\n";
            #my $start_diff = abs($transcript->seq_region_start() - $mane_transcript->seq_region_start());
            #print $start_diff, "\n";
            #my $end_diff = abs($transcript->seq_region_end() - $mane_transcript->seq_region_end());
            #if ($start_diff > 10 || $end_diff > 10) {
            #}
                #print $gene->stable_id(), "\t", $transcript->stable_id(), "\t", $transcript->seq_region_name(), "\t", $transcript->seq_region_start(), "\t", $transcript->seq_region_end(), "\n";
        #}
    }
}


##Checking the common genes between the two dictionaries just to make sure that 
## The code is actually doing what I think it is doing 
# Initialize an empty hash
my %common_genes;

# Loop through all gene names (keys) in the mane_transcripts hash
foreach my $gene_id (keys %mane_transcripts) {
    # Check if the gene is also in the non_mane_transcripts hash
    if (exists $non_mane_transcripts{$gene_id}) {
        # Store the gene and its transcripts in the common_genes hash
        $common_genes{$gene_id} = {
            'MANE_Transcripts' => $mane_transcripts{$gene_id},
            'Non_MANE_Transcripts' => $non_mane_transcripts{$gene_id}
        };
    }
}

my $diff_limit = 10;  # Define the default difference limit
GetOptions('diff_limit=i' => \$diff_limit) or die "Usage: $0 --diff_limit <number>\n"; #set diff limit, default = 10

open my $start_diff_file, '>', "start_diff_${diff_limit}bp.txt" or die "Could not open start_diff.txt: $!";
open my $end_diff_file, '>', "end_diff_${diff_limit}bp.txt" or die "Could not open end_diff.txt: $!";

# Initialize a hash to store the count of transcripts with an overhang for each gene
my %gene_transcript_counts;

# Loop through all gene names (keys) in the common_genes hash
foreach my $gene_id (keys %common_genes) {
    # Get the MANE and non-MANE transcript IDs for the gene
    my $mane_transcript_ids = $common_genes{$gene_id}{'MANE_Transcripts'};
    my $non_mane_transcript_ids = $common_genes{$gene_id}{'Non_MANE_Transcripts'};

    # Compare each MANE transcript to each non-MANE transcript
    foreach my $mane_transcript_id (@{$mane_transcript_ids}) {
        foreach my $non_mane_transcript_id (@{$non_mane_transcript_ids}) {
            # Fetch the transcript objects from the Ensembl database
            my $mane_transcript = $transcript_adaptor->fetch_by_stable_id($mane_transcript_id);
            my $non_mane_transcript = $transcript_adaptor->fetch_by_stable_id($non_mane_transcript_id);

            # Determine the strand of the transcript
            my $strand = $mane_transcript->strand();

            # Calculate the differences in the start and end sites
            my $start_diff;
            my $end_diff;
            if ($strand == 1) {
                $start_diff = $mane_transcript->seq_region_start() - $non_mane_transcript->seq_region_start();
                $end_diff = $mane_transcript->seq_region_end() - $non_mane_transcript->seq_region_end();
                # Count the number of transcripts that are 5' or 3' of the mane start site
                if ($start_diff > 0) {
                    if ($start_diff < $diff_limit) {
                        $count_5_prime++;
                        $gene_transcript_counts{$gene_id}++;
                        print $start_diff_file "$gene_id\t" . $non_mane_transcript->stable_id() . "\t" . $non_mane_transcript->seq_region_name() . "\t" . $non_mane_transcript->seq_region_start() . "\t" . $non_mane_transcript->seq_region_end() . "\t" . $non_mane_transcript->external_name() . "\n";
                    }
                } 
                elsif ($start_diff < 0) {
                    if (abs($start_diff) < $diff_limit) {
                        $count_3_prime++;
                        print $start_diff_file "$gene_id\t" . $non_mane_transcript->stable_id() . "\t" . $non_mane_transcript->seq_region_name() . "\t" . $non_mane_transcript->seq_region_start() . "\t" . $non_mane_transcript->seq_region_end() . "\t" . $non_mane_transcript->external_name() . "\n";
                    }
                }
            } else {
                $start_diff = $mane_transcript->seq_region_end() - $non_mane_transcript->seq_region_end();
                $end_diff = $mane_transcript->seq_region_start() - $non_mane_transcript->seq_region_start();
                # Count the number of transcripts that are 5' or 3' of the mane start site
                if ($start_diff < 0) {
                    if (abs($start_diff) < $diff_limit) {
                        $count_5_prime++;
                        $gene_transcript_counts{gene_id}++;
                        print $start_diff_file "$gene_id\t" . $non_mane_transcript->stable_id() . "\t" . $non_mane_transcript->seq_region_name() . "\t" . $non_mane_transcript->seq_region_start() . "\t" . $non_mane_transcript->seq_region_end() . "\t" . $non_mane_transcript->external_name() . "\n";
                    }
                } 
                elsif ($start_diff > 0) {
                    if ($start_diff < $diff_limit) {
                        $count_3_prime++;
                        print $start_diff_file "$gene_id\t" . $non_mane_transcript->stable_id() . "\t" . $non_mane_transcript->seq_region_name() . "\t" . $non_mane_transcript->seq_region_start() . "\t" . $non_mane_transcript->seq_region_end() . "\t" . $non_mane_transcript->external_name() . "\n";
                    }
                }
            }
        }
    }
}


open(my $number, '>', "rankedUpstream_${diff_limit}.txt") or die "Could not open file 'output.txt' $!";
# Remove the header from the hash
delete $gene_transcript_counts{'gene_id'};
## Counts number of transcripts per gene and ranks by highest number of transcripts per gene 
foreach my $gene_id (sort { $gene_transcript_counts{$b} <=> $gene_transcript_counts{$a} } keys %gene_transcript_counts) {
    print $number "$gene_id has $gene_transcript_counts{$gene_id} transcripts with a 5' overhang\n";
}
close $number;

# Open the file in write mode
open(my $fh, '>', "numberOftranscripts5prime_${diff_limit}.txt") or die "Could not open file 'output.txt' $!";

# Print the data to the file
print $fh "Number of transcripts upstream of the 5' mane TSS: $count_5_prime\n";
print $fh "Number of transcripts downstream of the 5' mane TSS: $count_3_prime\n";

# Close the file
close $fh;


# Close the files
close $start_diff_file;
close $end_diff_file;

# Now the common_genes hash contains the genes and their transcripts that are present in both hashes


# Print the MANE Select transcripts
#print "MANE Select Transcripts:\n";
#print Dumper(\%mane_transcripts);

# Print a newline for readability
#print "\n";

# Print the non-MANE Select transcripts
#print "Non-MANE Select Transcripts:\n";
#print Dumper(\%non_mane_transcripts);

##Print the overlapping genes
#print "OVERLAPPING GENES \n";
#print Dumper(\%common_genes);

