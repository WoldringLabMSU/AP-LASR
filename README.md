
**James VanAntwerp
10-07-23**


**AP-LASR: Automated Protein Libraries from Ancestral Sequence Reconstruction**

This code is used to automate the curation of combinatorial protein libraries from an automated ancestral sequence reconstruction workflow. Given a protein of interest, AP-LASR will curate a dataset of modern homologs, conduct indel-aware ancestral sequence reconstruction, and generate combinatorial libraries from ancestral proteins. This is valuable as a means to generate a diverse library of highly functional, stable protein sequences in an automated, objective manner. See Figure 1 for a representation of the software’s goals.

![FIg1](https://github.com/jjvanantwerp/Automated-ASR/assets/73084333/ce2704c1-62f0-469b-8a95-d6298a0d266e)
**Figure 1: A representation of the goals and motivations of AP-LASR. A) Phylogenetic trees relate modern homologs, but the sequences at internal nodes can be statistically reconstructed. B) Prediction of ancestral sequences comes with statistical uncertainty. Multiple amino acids may exist at a given position in the ancestral sequence. C) Most sites in an ancestral sequence (node) have very high confidence in one amino acid, but there are certain low-confidence positions which look more like the representation in B. These sites are distributed throughout the sequence. D) Combinations of alternative ancestral states (using lower confidence amino acids at a particular position rather than the most likely amino acid) can be used to generate libraries of sequences, where the libraries vary in size depending on the confidence needed for an amino acid to be included as an alternative sequence (confidence threshold). Because the libraries are generated combinatorially, they grow very rapidly with a decreasing confidence threshold. Different sequences with different numbers of low-confidence sites, have libraries which grow faster or slower with confidence threshold.**

**Build Status:**

This project was initially built as a personal tool to speed up the evaluation of amenable proteins for in-depth ASR. As such, it is not actively maintained and bug reports may take a long time to get fixed as I am not an experienced programmer and maintain this code as a hobby of sorts. The interactions between the Python Script and the OS are currently setup for Linux systems, not Windows, as it is best to run AP-LASR on a computing cluster rather than on a personal computer.

Input of multiple sequences is a supported feature, and often useful, but it is also prone to more errors as this feature has had limited testing and use. Controlling final dataset size is especially difficult for multi-sequence inputs; if your datasets (Final_Sequences.fasta) are too large, try manually adjusting the number of times that the alignment is cleaned by using option “-MSI’ (stands for multi-sequence iterations) and specifying a number of iterations (default 2). Additional known issues include access to NCBI’s BlastP through the BioPython API tool. Too many API calls in a short amount of time (i.e. if multiple jobs are running in parallel making calls from the same IP address) may result in that IP address to be temporarily blocked by NCBI, causing runs to end in an error state. When many requests are submitted over a longer period of time, such as when the low-homology sequence supplementation feature is turned on, this may result in throttling of requests, significantly lengthening the amount of time for runs to complete.


**Code Style:**

This project has no code style. Sorry! It’s all written in my own messed up way. I did run it through PEP8 formatting at least. I’m not an animal.

**Features:**

AP-LASR can run in one of four modes - “Assist” mode gives more detailed help than just “-h” would. “ASR” mode is the majority of intended use - it will take in user sequences and other parameters, then run ASR and build libraries. See Figure 2 for a representation of the workflow of AP-LASR in ASR mode. “MakeFigures” mode will generate figures for a completed ASR run, given the name of that directory. If this feature is to be used, the python module Matplotlib is required. “RemakeLibraries” mode will regenerate the libraries from a completed ASR run with a user-specified threshold confidence, expressed as a decimal. 

![Fig2](https://github.com/jjvanantwerp/Automated-ASR/assets/73084333/64e714aa-b728-426e-8234-bd97bb2cb532)
**Figure 2: A representation of the workflow of AP-LASR (Low_Similarity Supplementation is an optional step)**


**Example Inputs:**

python3 APLASR.py ASR -i input.fasta -n ASR_RUN_10_7_23 -s 3000 --iqtree iqtree2 —cdhit CDHIT
	
The above line runs APLASR in ASR mode, pulling input sequences from the file ‘input.fasta’, and stores the output in a directory called “ASR_RUN_10_7_23”. It caps the final dataset to 3,000 sequences, and is running in an environment where “iqtree2” is the executable for IQ-Tree and “CDHIT” is the executable for CD-Hit.

 
python3 APLASR.py RemakeLibraries -n ASR_RUN_10_7_23 -t 0.085

The above line runs APLASR in RemakeLibraries mode, regenerating libraries for the above reconstruction with a confidence threshold of 8.5%.

 
python3 APLASR.py MakeFigures -n ASR_RUN_10_7_23
 
The above line runs APLASR in MakeFigures mode, generating figures for the above reconstruction.


python3 APLASR.py ASR -i input.fasta

The above line runs APLASR in ASR mode, with the input sequence in the file “input.fasta”. All other parameters are default.

 
python3 APLASR.py ASR -i MRFPSIFTAVLFAASSALAAPVNTTTEDETAQIPAEAVIGYLDLEGDFDVAVLPFSNSTNNGLLFINTTIASIAAKEEGVSLDKREAEAWHWLQLKPGQPMYKREAEAEAWHWLQLKPGQPMYKREADAEAWHWLQLKPGQPMYKREADAEAWHWLQLKPGQPMY -n alpha_ASR
  
The above line runs APLASR in ASR mode, with the input sequence as the list of amino acids pasted directly into the command line (in this case, yeast mating factor alpha-1 from S. Cerevisiae, UniProt P01149). The directory name is specified as “alpha_ASR” and all other parameters are left default.


**How to Use:**

It is recommended that AP-LASR is run on a computing cluster, as it is computationally intensive, and can run for dozens of hours. The set-up descriptions assume that is the case for your use as well. Download the python file in this GitHub repository, and upload or copy it into the computing cluster. In the same directory, create a file to contain the input sequence(s), in fasta format. Determine how you can get access to the three software dependences (IQ-Tree, CD-Hit, MAFFT) on the computing cluster - cluster administrators should be able to help with this, if there is a built-in module or of they need to add one for you. Additionally, ensure that you have access to the Biopython module. Once the proper modules are loaded, AP-LASR can be run from the command line as described in the “Example Inputs” section of this README. The time required to run AP-LASR depends significantly on the length of and number of sequences in the final sequence alignment. Sequence supplementation can also add significantly to the runtime, depending on the number of sequences which are resubmitted. Once AP-LASR has finished, you can get a description of the file locations and contents by running AP-LASR in “Assist” mode. If you wish to use AP-LASR locally, then go and download these three software: IQ-Tree, CD-Hit, and MAFFT (links below, functional as of 10-07-23). Add their executables to your PATH variable, and run AP-LASR from the command line the same as above. 

CD-Hit:  https://github.com/weizhongli/cdhit/releases 

MAFFT:  https://mafft.cbrc.jp/alignment/software/

IQ-Tree:  http://www.iqtree.org/#download


**Contributors:**

James VanAntwerp wrote the lions’ share of the code, and did all the bug-squashing, documentation, and testing, and wrote the manuscript drafts.
Daniel Woldring was research advisor, idea-generator, and project manager. He also provided significantly to the manuscript in terms of content and editing.
Ben Dolgikh and Patrick Finneran contributed some code for specific steps of the workflow, and strategic advice.
Mehrsa Mardikoraem contributed statistical evaluations of test conditions, as well as contributions to and editing of the manuscript and figures.
Nathan Pascual made contributions to and editing of the manuscript and figures.
Many thanks to all of them for their support in this project, and to my freinds and family for encouraging me through the long and sometimes difficult process that was the creation of AP-LASR.


**License:**

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
