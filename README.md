# BLAST_Genes
This python script BLASTs standard gene sequences against a set of genomes to determine which genomes have the genes.

Files: blast_genes.py: python script 

Installation

Install Python Download from python.org. Make sure to check "Add Python to PATH" during installation.

Install Required Python Packages biopython pandas click tqdm

Install NCBI BLAST+ Download and install from: https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ Add the BLAST+ folder to your system PATH.

Setup for Command Line Access (Optional but Recommended) On Windows Make sure the folder where your blast_genes.py file is located is in your system PATH. To find the path, open the folder and copy the address bar (e.g. C:\Users\YourName\scripts), then: Press ⊞ Win, search for "Environment Variables" Under "System Variables", find Path, click Edit → New → paste the folder path. Click OK and restart your terminal.

On Linux If using pip install --user, ensure this is in your shell config (~/.bashrc, ~/.zshrc, etc.): export PATH="$HOME/.local/bin:$PATH" Run: source ~/.bashrc # or ~/.zshrc

Usage python blast_genes.py -f genome_sequences.fasta -d Standard_gene_sequences.fasta 
All input files must be in FASTA format. BLAST databases are built automatically.

Output Files: FinalResults.csv 
All successful gene matches at least 500 bp long.

Genes must bind to sequence with a binding length of ≥500 bp and match ≥88% identity. Output folders are overwritten each run — back up results if needed.
