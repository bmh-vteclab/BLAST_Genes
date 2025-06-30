import os
import platform
import subprocess
import click
import pandas as pd
from Bio import SeqIO
from tqdm import tqdm

# Minimum allowed alignment length
MIN_GENE_LENGTH = 500

@click.command(help="blast_genes: BLASTs standard gene sequences against a set of genomes to determine which genomes have the genes.")
@click.option('-f', '--fasta', type=click.Path(exists=True), help='Path to the sequences of genomes to be examined fasta file.')
@click.option('-d', '--gene', type=click.Path(exists=True), help='Path to the standard gene sequences fasta file.')

def main(fasta, gene):
    system = platform.system().lower()
    check_blast_installed(system)

    sequences = list(SeqIO.parse(fasta, 'fasta'))

    make_blast_db(gene)
    results = []

    for seq_record in tqdm(sequences, desc="Processing sequences"):
        seq_id = seq_record.id
        temp_query = 'temp_query.fasta'
        with open(temp_query, 'w') as f:
            SeqIO.write(seq_record, f, 'fasta')

       gene_hits = blast_query(system, temp_query, 'gene_db')

        #matched_result = []
        
        for gene in gene_hits.itertuples(index=False):
            row = [seq_id, gene.sacc, gene.length, gene.qstart, gene.qend, gene.pident]
            results.append(row)
  
    os.remove(temp_query)
    
    columns = ['Strain', 'Gene', 'gene length', 'gene Start', 'gene End', 'gene %indent']

    results_df = pd.DataFrame(results, columns=columns)
    results_df.to_csv('FinalResults.csv', index=False)
    
    cleanup_temp_files()
    click.echo("\nAll finished!")

def check_blast_installed(system):
    try:
        command = ['blastn', '-version'] if system != 'windows' else 'blastn -version'
        subprocess.run(command, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True, shell=(system=='windows'))
    except (subprocess.CalledProcessError, FileNotFoundError):
        click.echo("BLAST+ is not installed or not added to PATH.")
        raise SystemExit(1)

def make_blast_db(gene_file):
    subprocess.check_output(f'makeblastdb -in "{gene_file}" -parse_seqids -dbtype nucl -out gene_db -title gene_db', shell=True)

def blast_query(system, query_file, db_name):
    output_file = 'temp_result.csv'
    if system == 'windows':
        subprocess.check_output(f'blastn -task blastn -query "{query_file}" -db "{db_name}" -out "{output_file}" -outfmt "6 sacc qseqid qstart qend length pident sstrand" -word_size 7 -evalue 1000 -perc_identity 88', shell=True)
    else:
        subprocess.check_output(f'blastn -task blastn -query {query_file} -db {db_name} -out {output_file} -outfmt "6 sacc qseqid qstart qend length pident sstrand" -word_size 7 -evalue 1000 -perc_identity 88', shell=True)

    if os.path.exists(output_file):
        df = pd.read_csv(output_file, sep='\t', header=None, names=['sacc', 'qseqid', 'qstart', 'qend', 'length', 'pident', 'sstrand'])
        os.remove(output_file)
        return df[df.length >= MIN_GENE_LENGTH]
    else:
        return pd.DataFrame(columns=['sacc', 'qseqid', 'qstart', 'qend', 'length', 'pident', 'sstrand'])

def cleanup_temp_files():
    for prefix in ['gene_db']:
        for ext in ['nhr', 'nin', 'nsq', 'nos', 'njs', 'not', 'ntf', 'nto', 'ndb', 'nog']:
            file = f'{prefix}.{ext}'
            if os.path.exists(file):
                os.remove(file)

if __name__ == '__main__':
    main()