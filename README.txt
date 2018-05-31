Copyright Etienne Muller (2016)

biomol.cfb@gmail.com

outLyzer is a computer program whose purpose is to detect variations,
specifically low allele frequency variation, in next generation 
sequencing data (tumor samples, mosaïc mutation)

This software is a system to highlight mutations that must be used 
with caution. We can not guarantee the accuracy of informations and
predictions provided.

This software is governed by the CeCILL license under French law and
abiding by the rules of distribution of free software.  You can  use, 
modify and/ or redistribute the software under the terms of the CeCILL
license as circulated by CEA, CNRS and INRIA at the following URL
"http://www.cecill.info". 

As a counterpart to the access to the source code and  rights to copy,
modify and redistribute granted by the license, users are provided only
with a limited warranty  and the software's author,  the holder of the
economic rights,  and the successive licensors  have only  limited
liability. 

In this respect, the user's attention is drawn to the risks associated
with loading,  using,  modifying and/or developing or reproducing the
software by the user in light of its specific status of free software,
that may mean  that it is complicated to manipulate,  and  that  also
therefore means  that it is reserved for developers  and  experienced
professionals having in-depth computer knowledge. Users are therefore
encouraged to load and test the software's suitability as regards their
requirements in conditions enabling the security of their systems and/or 
data to be ensured and,  more generally, to use and operate it in the 
same conditions as regards security. 

The fact that you are presently reading this means that you have had
knowledge of the CeCILL license and that you accept its terms.

outLyzer version 2.0 31/05/2018

Author: Etienne Muller
E-mail: biomol.cfb@gmail.com
Sources: http://github.com/EtieM/outLyzer



#               _     __                    
#    ___  _   _| |_  / / _   _ _______ _ __ 
#   / _ \| | | | __|/ / | | | |_  / _ \ '__|
#  | (_) | |_| | |_/ /__| |_| |/ /  __/ |   
#   \___/ \__,_|\__\____/\__, /___\___|_|   
#                        |___/           



CONTENTS OF THIS FILE
---------------------

* Introduction
* Requirements
* Installation
* Utilisation




INTRODUCTION
------------

OutLyzer is a variant-caller conceived for low allele-ratio mutations
detection, based on sequencing background noise evaluation. It evaluates if
the mutation is significantly different from background noise, using modified
Thompson tau technique.
This program was conceived in the Department of Cancer Biology and Genetics,
Francois Baclesse Cancer Center, Caen, France 
Program can be downloaded at: http://github.com/EtieM/outLyzer


REQUIREMENTS
------------
- Linux OS
- Python v2.7
- Python librairies: subprocess / numpy / scipy / argparse / multiprocessing
- Samtools v1.2 or v1.3


INSTALLATION
------------

Before launching outLyzer program, it is preferable to set Samtools path as an environment variable:
- For one session: $export samtools=/path/to/samtools
- Permanantly: write previous command line in .bashrc file

Otherwise complete Samtools path must be indicated in outLyzer launching command line.



UTILISATION
-----------
outLyzer has 2 main operating modes:

- calling mode: $python outLyzer.py calling -bed regionFile.bed -ref referenceFile.fa -bam fileToAnalyse.bam -output /path/to/resultDir/[-arguments]
				This mode works as a standard variant-caller and analyse a whole BAM file to highlight mutations, compiled in a VCF file.

				optional arguments:
				  -h, --help            show this help message and exit
				  -samtools SAMTOOLS    Complete Samtools path if not specified in environment
										variable
				  -pythonPath PYTHONPATH
										Complete python path if different from default python
										version
				  -core CORE            define number of cores used to process analysis [1]
				  -cut CUT              defines into how many parts bed file is divide [3]
				  -bed BED              bed file required for analysis [REQUIRED]
				  -bam BAM              bam File to analyze [REQUIRED]
				  -ref REF              faidx indexed reference sequence file (fasta)
										[REQUIRED]
				  -output OUTPUT        output Path To write results [REQUIRED]
				  -t T                  Student t value used in modified Thompson tau
										technique [0.001]
				  -bal BAL              minimum Forward / Reverse read proportion [0.3]
				  -Q Q                  minimum average Phred Score to be considered as a real
										mutation (only relevant for SNP) [20]
				  -SDQ SDQ              maximum Standard deviation authorized for average
										Phred Score [7]
				  -WS WS                Window Size: region (number of bp) around the mutation
										on which background noise have to be determined [200]
				  -WSmin WSMIN          Window Size Minimum Size: minimum region size (number
										of bp) required for analysis [10]
				  -x X                  Multiplicative factor that specifies how often the
										mutation must be above background noise [2]
				  -AS                   Analysis sensitivity: Returns an additional file
										containing analysis average sensitivity for each line
										of bed file
				  -FRcor                Forward Reverse Correction: take into account any
										imbalance in the Forward-Reverse reads distribution in
										the Forward / Reverse alternative Read Proportion
										(-bal option)
				  -HSM HSM              HotSpot Metrics: Produce sensitivity Threshold for
										HotSpot positions, in an additional file. Requires
										formated HotSpot File in argument (see documentation
										for more details).
				  -verbose VERBOSE      If verbose mode is set to 1, details analysis process
										steps [0]

				
				Precisions for HSM option:
				/!\ File required for HotSpot Metrics must be formatted as follows:
				
				chrN	startPosition	Annotation
				
				Each column must be separated by a tabulation, and annotation column must be present.
				ex: chr12	25398284	KRAS_codon12
				
				
				It will return a tabulated file containing for each position a local estimation of sensitivity, displayed as a percentage.

				
-positionAnalysis mode: $ python outLyzer.py positionAnalysis -bam fileToAnalyse.bam -ref referenceFile.fa -position chr12:123456 [-arguments]
						
						This mode gives an evaluation of sequencing data and local noise background for one chromosomic position
						
						
						optional arguments:
						  -h, --help          show this help message and exit
						  -samtools SAMTOOLS  Complete Samtools path if not specified in environment
											  variable
						  -bam BAM            bam File to analyze [REQUIRED]
						  -position POSITION  chromosomic position to analyze (ex: chr3:123456789)
											  [REQUIRED]
						  -ref REF            faidx indexed reference sequence file (fasta) [REQUIRED]
						  -t T                Student t value used in modified Thompson tau technique
											  [0.001]
						  -bal BAL            minimum Forward / Reverse read proportion [0.3]
						  -Q Q                minimum average Phred Score to be considered as a real
											  mutation (only relevant for SNP) [20]
						  -SDQ SDQ            maximum Standard deviation authorized for average Phred
											  Score [7]
						  -WS WS              Window Size: region (number of bp) around the mutation
											  on which background noise have to be determined [200]
											  
						It directly displays results in standard output as follows:
						
						['1', '0', '0', '0', '2', '0', '1', '0', '0', '0', '1', '1', '0', '1',
						'1', '0', '0', '1', '0', '0', '0', '1', '0', '2', '1', '1', '0', '1', 
						'0', '0', '2', '2', '0', '0', '0', '1', '0', '1', '0', '0', '0', '0', 
						'0', '0', '0', '0', '1', '2', '0', '0', '1', '1', '0', '0', '1', '0', 
						'2', '0', '0', '1', '0', '1', '0', '0', '0', '1', '0', '3', '0', '0', 
						'0', '1', '0', '2', '0', '0', '0', '2', '0', '0', '0', '1', '2', '0', 
						'0', '0', '0', '2', '1', '2', '1', '0', '0', '3', '0', '0', '0', '0', 
						'0', '1', '964', '1', '0', '0', '0', '0', '0', '0', '0', '1', '0', '1',
						'1', '0', '0', '0', '1', '2', '0', '0', '0', '0', '0', '1', '5', '1', '0',
						'1', '0', '0', '0', '0', '1', '0', '1', '0', '1', '1', '0', '0', '0', '2', 
						'0', '0', '0', '0', '2', '1', '1', '5', '0', '1', '2', '2', '0', '0', '1', 
						'0', '1', '0', '0', '1', '0', '1', '1', '4', '0', '0', '0', '2', '0', '0', 
						'1', '0', '0', '2', '2', '1', '2', '2', '1', '3', '0', '5', '4', '1', '0', 
						'1', '2', '0', '0', '2', '0', '1', '0', '2', '0', '0', '0', '4', '1']

						Mutation Position:         chr12:25378562
						Reference Allele:          C
						Alternative Allele:        T
						Depth:                     1835
						Allele Frequency (%):     52.48
						Phred Quality:             26.5
						Phred Standard Deviation:  1.73
						Forward / Reverse alt:     494/469     (51.3% / 48.7%)
						overAll Balance:           50.5% / 49.5% 
						Corrected alt F/R:         51.0% / 49.0% 
						Raw background Noise:      3
						Stretch nearby:            None
						Motif nearby:              None

						
						The sequence of numbers represents, for each genomic position in the analyzed window,
						the number of alternative reads.
						
						Alternative Allele: mutated base. (WT = Wild Type - No mutation on this position)
						Depth = total number of reads aligned to this position
						Allele Frequency = proportion of alternative reads
						Phred Quality = average PHRED quality for all alternative reads
						Phred standard deviation = standard deviation for average Phred quality score mentioned above
						Forward / Reverse = number of alternative reads sequenced in forward / reverse
						overAll Balance: proportion of reads (wild-type AND mutated) sequenced in forward / reverse
						Corrected alt F/R: ponderation of alternative Forward / Reverse balance based on overall balance
						Raw background = sequencing background noise around the mutation, expressed in number of reads
						Stretch nearby: indicates if there is a strecth nearby the mutation
						Motif nearby: indicates if there is a repetitive DNA-sequence motif nearby the mutation
						
						
-Utilisation test:

$ python outLyzer_v2.0.py -ref reference.fasta -bed test_dataSet.bed -bam  test_dataSet.bam -HSM HotSpot_positions_HSM.bed -AS -FRcor -output outputPath

Results should correspond to resultsExample files
