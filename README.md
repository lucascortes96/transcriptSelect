# transcriptSelect
Will become part of a GENCODE project that works to pick transcripts that can be extended or clipped at the 5' TSS or the 3' transcription end automatically using an algorithm. Will eventually pull in CAGE, RNA Seq and other data types in order to support reasons for the changes. 


## retrieve5prime.pl
Grabs transcripts from protein coding genes that are within 'n' base pairs of the 5' TSS of the MANE select transcript. Specifically 'rankedUpstream.txt' provides the top genes that have the highest number of transcripts upstream of the 5' MANE TSS that fall within the 'n' base pair limit. 

##parseOutput.py
Simply grabs a couple of general stats from the retrieve5prime script. 
