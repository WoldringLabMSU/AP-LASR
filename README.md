# AP-LASR: Automated Protein Libraries from Ancestral Sequence Reconstruction

AP-LASR automates the curation of combinatorial protein libraries from an automated ancestral sequence reconstruction workflow. Given a protein of interest, it will curate a dataset of modern homologs, conduct indel-aware ancestral sequence reconstruction, and generate combinatorial libraries from ancestral proteins. This is valuable as a means to generate a diverse library of highly functional, stable protein sequences in an automated, objective manner.

For more information about the motivation behind AP-LASR, check out the pre-print: https://www.biorxiv.org/content/10.1101/2023.10.09.561537v1

>![FIg1](https://github.com/jjvanantwerp/Automated-ASR/assets/73084333/ce2704c1-62f0-469b-8a95-d6298a0d266e)
>**Figure 1:** A representation of the goals and motivations of AP-LASR. **A)** Phylogenetic trees relate modern homologs, but the sequences at internal nodes can be statistically reconstructed. **B)** Prediction of ancestral sequences comes with statistical uncertainty. Multiple amino acids may exist at a given position in the ancestral sequence. **C)** Most sites in an ancestral sequence (node) have very high confidence in one amino acid, but there are certain low-confidence positions which look more like the representation in B. These sites are distributed throughout the sequence. **D)** Combinations of alternative ancestral states (using lower confidence amino acids at a particular position rather than the most likely amino acid) can be used to generate libraries of sequences, where the libraries vary in size depending on the confidence needed for an amino acid to be included as an alternative sequence (confidence threshold). Because the libraries are generated combinatorially, they grow very rapidly with a decreasing confidence threshold. Different sequences with different numbers of low-confidence sites, have libraries which grow faster or slower with confidence threshold.


# Features

AP-LASR can run in one of four modes: 
* `ASR` mode is the majority of intended use - it will take in user sequences and other parameters, then run ASR and build libraries. See Figure 2 for a representation of the workflow of AP-LASR in ASR mode.
+ `Assist` mode gives more detailed help than just `-h` would. 
* `MakeFigures` mode will generate figures for a completed ASR run, given the name of that directory. If this feature is to be used, the python module Matplotlib is required. 
+ `RemakeLibraries` mode will regenerate the libraries from a completed ASR run with a user-specified threshold confidence, expressed as a decimal. 

>![Fig2](https://github.com/jjvanantwerp/Automated-ASR/assets/73084333/64e714aa-b728-426e-8234-bd97bb2cb532)
**Figure 2:** Workflow of AP-LASR (Low_Similarity Supplementation is an optional step).

# Using AP-LASR

## Installing AP-LASR

1. Download the python file in this GitHub repository.
2. In the same directory, create a file to contain the input sequence(s), in fasta format. 
3. Install the following dependencies: IQ-Tree, CD-Hit, MAFFT, and Biopython. If you do not have install priviliges, contact your compute cluster to help with this. Built-in modules or virtual environments may be helpful in these cases. 

    * **CD-Hit:**  https://github.com/weizhongli/cdhit/releases 

    * **MAFFT:**  https://mafft.cbrc.jp/alignment/software/

    * **IQ-Tree:**  http://www.iqtree.org/#download

    While not necessary to run AP-LASR, Matplotlib should be installed to generate figures.

    * **(Optional)Matplotlib:** https://matplotlib.org/stable/users/getting_started/ 

4. Once the proper modules are loaded, AP-LASR can be run from the command line as described in the following section of this README. 


## Basic Usage

`python3 APLASR.py ASR -i input.fasta -n ASR_RUN_10_7_23 -s 3000 --iqtree iqtree2 —cdhit CDHIT`
	
The above line runs APLASR in ASR mode, pulling input sequences from the file ‘input.fasta’, and stores the output in a directory called “ASR_RUN_10_7_23”. It caps the final dataset to 3,000 sequences, and is running in an environment where “iqtree2” is the executable for IQ-Tree and “CDHIT” is the executable for CD-Hit.

## Options/Flags
For the default `ASR` mode:

`-i` 
- specifies the name of the input Fasta file, or can take a raw protein sequence. This is the only mandatory input

`-n` 
- specifies the name of the output directory. The default is "ASR".

`-s` 
- specifies the desired maximum size of the final dataset of modern homologs. This allows users some control over the time AP-LASR will take to run, and the detail/quality of the ASR. The default is 500.

`-supplement` 
- will turn on supplementation at the provided similarity cutoff (entered as a decimal). If you do turn on supplementation, there is not a default value but we reccomend users specify something between 0.7 and 0.85. 

`-iqtree` 
- allows users to specify an executable for IQTree other than the default ("iqtree").
    
`-cdhit`
- allows users to specify an executable for CH-Hit other than the default ("cd-hit").

`-mafft` 
- allows users to specify an executable for MAFFT other than the default ("mafft").

`-MSI` 
- allows users to specify the number of times the function "Post_MAFFT_Processing" iterates over the raw alignment to generate the final dataset. Each iteation will remove a large number of sequences, more if the different input sequecnes are dissimilar. Default is 2.

`-t` 
-   specifies the library generation threshold, expressed as a decimal. 

`-o`, or `--outgroup`
-  allows users to provide an outgroup via a FASTA file, which will be included in the final sequence alignment for tree rooting, etc.

As briefly described above, there are three other modes for AP-LASR, with each having specific options.

`RemakeLibraries`
- This mode will generate new combinatorial libraries from ancestral proteins of a previous ASR run with a specified threshold confidence. This will allow generation of a library with a different size than what was generated with the default threshold confidences. Specify the directory where the ASR results of interest are stored.
-  `RemakeLibraries`-specific options:
    - specify the directory where the results of interest are stored with the `-n` flag.
    - a new threshold may be defined with the `-t` option
    - `-A` or `--AltAll` 
        -  will make libraries in an "Alt-All" mode, which gives one highly mutated alternative sequence for each ancestor. No threshold should be specified if this mode is used.

`MakeFigures` 
- This mode will generate figures from ASR data of a previous ASR run
- specify the directory where the results of interest are stored with the `-n` flag.

`Assist` 
- This mode gives more detailed help than just `-h` would. 

## Helpful Examples
 
`python3 APLASR.py RemakeLibraries -n ASR_RUN_10_7_23 -t 0.085`

The above line runs APLASR in RemakeLibraries mode, regenerating libraries for the above reconstruction with a confidence threshold of 8.5%.

 
`python3 APLASR.py MakeFigures -n ASR_RUN_10_7_23`
 
The above line runs APLASR in MakeFigures mode, generating figures for the above reconstruction.

`python3 APLASR.py ASR -i input.fasta`

The above line runs APLASR in ASR mode, with the input sequence in the file “input.fasta”. All other parameters are default.

`python3 APLASR.py ASR -i MRFPSIFTAVLFAASSALAAPVNTTTEDETAQIPAEAVIGYLDLEGDFDVAVLPFSNSTNNGLLFINTTIASIAAKEEGVSLDKREAEAWHWLQLKPGQPMYKREAEAEAWHWLQLKPGQPMYKREADAEAWHWLQLKPGQPMYKREADAEAWHWLQLKPGQPMY -n alpha_ASR`
  
The above line runs APLASR in ASR mode, with the input sequence as the list of amino acids pasted directly into the command line (in this case, yeast mating factor alpha-1 from S. Cerevisiae, UniProt P01149). The directory name is specified as “alpha_ASR” and all other parameters are left default.

`python3 APLASR.py ASR -i input.fasta -o outgroup.fasta`

The above line runs APLASR in ASR mode, with the input sequence in the file “input.fasta” and an outgroup sepcified in "outgroup.fasta". All other parameters are default.

# Helpful Tips

### Using Computing Clusters is Recommended

It is recommended that AP-LASR is run on a computing cluster, as it is computationally intensive, and can run for dozens of hours. The set-up descriptions assume that is the case for your use as well.

If you wish to use AP-LASR locally, ensure the three software: IQ-Tree, CD-Hit, and MAFFT. Add their executables to your PATH variable, and run AP-LASR from the command line the same as above. 

### AP-LASR depends on sequence length and number of sequences aligned. 

The time required to run AP-LASR depends significantly on the length of and number of sequences in the final sequence alignment. Sequence supplementation can also add significantly to the runtime, depending on the number of sequences which are resubmitted. 

### Run AP-LASR in `Assist` Mode for a Description of Outputs 
Once AP-LASR has finished, you can get a description of the file locations and contents by running AP-LASR in `Assist` mode. 

### Code Style:

This project has no code style. Sorry! It’s all written in my own messed up way. I did run it through PEP8 formatting at least. I’m not an animal.

### Build Status:

This project was initially built as a personal tool to speed up the evaluation of amenable proteins for in-depth ASR. As such, it is not actively maintained and bug reports may take a long time to get fixed as I am not an experienced programmer and maintain this code as a hobby of sorts. The interactions between the Python Script and the OS are currently setup for Linux systems, not Windows, as it is best to run AP-LASR on a computing cluster rather than on a personal computer.

### Use multiple sequences with caution.

Input of multiple sequences is a supported feature, and often useful, but it is also prone to more errors as this feature has had limited testing and use. Controlling final dataset size is especially difficult for multi-sequence inputs; if your datasets (Final_Sequences.fasta) are too large, try manually adjusting the number of times that the alignment is cleaned by using option “-MSI’ (stands for multi-sequence iterations) and specifying a number of iterations (default 2). 

### Running multiple parallel AP-LASR jobs, errors may result due to NCBI API limits.
Additional known issues include access to NCBI’s BlastP through the BioPython API tool. Too many API calls in a short amount of time (i.e. if multiple jobs are running in parallel making calls from the same IP address) may result in that IP address to be temporarily blocked by NCBI, causing runs to end in an error state. When many requests are submitted over a longer period of time, such as when the low-homology sequence supplementation feature is turned on, this may result in throttling of requests, significantly lengthening the amount of time for runs to complete.

# Contributors

James VanAntwerp wrote the lions’ share of the code, and did all the bug-squashing, documentation, and testing, and wrote the manuscript drafts.

Daniel Woldring was research advisor, idea-generator, and project manager. He also provided significantly to the manuscript in terms of content and editing.

Ben Dolgikh and Patrick Finneran contributed some code for AP-LASR and provided strategic advice.

Mehrsa Mardikoraem contributed statistical evaluations of test conditions, as well as contributions to and editing of the manuscript and figures.

Nathan Pascual made contributions to and editing of the manuscript and figures.

Many thanks to all of them for their support in this project, and to my freinds and family for encouraging me through the long and sometimes difficult process that was the creation of AP-LASR.

# License

MIT License

Copyright (c) 2023 James Van Antwerp

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

# Change Log

James VanAntwerp, 10-07-23, page semi-finalized.

Nathan Pascual, 10-26-23, minor formatting edits and updates for Outgroup Functionality and "Alt-All" mode.