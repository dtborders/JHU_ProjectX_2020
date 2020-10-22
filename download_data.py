import os
import argparse
import numpy as np
import pandas as pd
from Bio import SeqIO, Entrez
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

'''
Usage:
    python download_data.py --save_dir {directory to save data, default is data}
        --input file {CSV file with all the host/phage id numbers}
'''

def main(save_dir, input_file, num_phages, num_hosts, email):
    Entrez.email = email

    comb_dsets = pd.read_csv(input_file)
    phages = list(set(comb_dsets['phage']))[1:]
    bacterias = list(set(comb_dsets['host']))[1:]

    # Download the Phages
    if num_phages is None or num_phages>=len(phages):
        handle = Entrez.efetch(db="nuccore", id=phages, rettype="gb", retmode="xml")
        phages_ass = phages
    else:
        handle = Entrez.efetch(db="nuccore", id=phages[0:num_phages], rettype="gb", retmode="xml")
        phages_ass = phages[0:num_phages]

    print("Downloading {} phages".format(len(phages_ass)))

    phages_download = Entrez.read(handle)
    handle.close()

    if not os.path.isdir(save_dir):
        os.mkdir(save_dir)

    phage_save_dir = os.path.join(save_dir, 'phages/')
    if not os.path.isdir(phage_save_dir):
        os.mkdir(phage_save_dir)



    for phage, phage_id in zip(phages_download, phages_ass):
        phage_id = str(phage_id)
        if not os.path.isdir(os.path.join(phage_save_dir, phage_id)):
            os.mkdir(os.path.join(phage_save_dir, phage_id))

        genes = np.array([foo for foo in phage['GBSeq_feature-table'] if foo['GBFeature_key'] == 'gene'])
        whole_seq = phage['GBSeq_sequence']

        gene_numbers = []
        for gene in genes:
            g_info = gene['GBFeature_intervals'][0]
            gene_numbers.append([int(g_info['GBInterval_from']), int(g_info['GBInterval_to'])])

        sequences = SeqRecord(
            Seq(whole_seq),
            id=phage_id,
            description="phage dna",
        )
        with open(os.path.join(phage_save_dir, str(phage_id), "{}.fasta".format(phage_id)), "w") as output_handle:
            SeqIO.write(sequences, output_handle, "fasta")

        np.save(os.path.join(phage_save_dir, str(phage_id), "{}.npy".format(phage_id)), np.array(gene_numbers))

    # Download the Hosts
    print("Downloading Hosts")
    host_save_dir = os.path.join(save_dir, 'hosts/')
    if not os.path.isdir(host_save_dir):
        os.mkdir(host_save_dir)

    bac_to_download = bacterias if (num_hosts is None or num_hosts>=len(bacterias)) else bacterias[0:num_hosts]
    handle = Entrez.efetch(db="nuccore", id=bac_to_download, rettype="gb", retmode="xml")
    host_infos = Entrez.read(handle)
    handle.close()

    for host_id, host_info in zip(bac_to_download, host_infos):

        host_id = host_id.split('.')[0]
        csvdir = os.path.join(host_save_dir, host_id)


        #two cases to deal with, sometimes nuccore has all the data needed, others it doesnt

        if 'GBSeq_sequence' in host_info.keys():
            print("\tDownloading {} with path 1".format(host_id))
            genes = np.array([foo for foo in host_info['GBSeq_feature-table'] if foo['GBFeature_key'] == 'gene'])
            whole_seq = host_info['GBSeq_sequence']
            gene_numbers = []
            for gene in genes:
                g_info = gene['GBFeature_intervals'][0]
                gene_numbers.append([int(g_info['GBInterval_from']), int(g_info['GBInterval_to'])])

            sequences = SeqRecord(
                Seq(whole_seq),
                id=host_id,
                description="host dna",
            )

            if not os.path.isdir(csvdir):
                os.mkdir(csvdir)

            with open(os.path.join(csvdir, "{}.fasta".format(phage_id)), "w") as output_handle:
                SeqIO.write(sequences, output_handle, "fasta")

            np.save(os.path.join(csvdir, "{}.npy".format(phage_id)), np.array(gene_numbers))
        else:
            print("\tDownloading {} with path 2".format(host_id))

            handle = Entrez.efetch(db="nuccore", id=host_id, rettype="fasta", retmode="xml")
            sequence = Entrez.read(handle)
            handle.close()

            whole_seq = sequence[0]['TSeq_sequence']

            handle2 = Entrez.esearch(db="gene", term=host_id, retmax=100000000, retmode="xml")
            genes = Entrez.read(handle2)
            handle2.close()

            if len(genes['IdList']) == 0:
                print("\t\tError, could not find any genes on ncbi gene database associates w/ host id ", host_id)
                continue

            handle3 = Entrez.efetch(db="gene", id=genes['IdList'], rettype='gene_table', retmode="xml")
            gene_infos = Entrez.read(handle3, validate=False)
            handle3.close()

            gene_locs = []
            for ii, cres in enumerate(gene_infos):
                #if host_id in cres['Entrezgene_gene-source']['Gene-source']['Gene-source_src-str1']:
                if len(cres['Entrezgene_locus'])<2:
                    #print('here')
                    continue
                seq_info = cres['Entrezgene_locus'][1]['Gene-commentary_seqs'][0]['Seq-loc_int']['Seq-interval']
                gene_locs.append([int(seq_info['Seq-interval_to']), int(seq_info['Seq-interval_from'])])
                #else:
                #    print(host_id, cres['Entrezgene_gene-source']['Gene-source']['Gene-source_src-str1'], 'invalid gene found')


            sequences = SeqRecord(
                Seq(whole_seq),
                id=host_id,
                description="host dna",
            )

            if not os.path.isdir(csvdir):
                os.mkdir(csvdir)

            with open(os.path.join(csvdir, "{}.fasta".format(host_id)), "w") as output_handle:
                SeqIO.write(sequences, output_handle, "fasta")

            np.save(os.path.join(csvdir, "{}.npy".format(host_id)), np.array(gene_locs))


    return 1

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--save_dir', default="data/", help="directory to store all the files")
    parser.add_argument('--input_file', default='Combined2.csv',
                        help="file to get phage/host id numbers from")
    parser.add_argument('--num_phages_to_download', default=None, type=int, help="the number of phages to downlaod, \
                                                                       if None will download all")
    parser.add_argument('--num_hosts_to_download', default=None, type=int, help="the number of hosts to downlaod, \
                                                                       if None will download all")
    parser.add_argument('--email', type=str, help="email address")
    args = parser.parse_args()
    main(args.save_dir, args.input_file, args.num_phages_to_download, args.num_hosts_to_download, args.email)
