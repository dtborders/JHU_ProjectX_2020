import argparse
import requests
import json
import os
import urllib.request
from Bio import SeqIO, Entrez
import os
import re
'''
Usage:
    python download_data.py --save_dir {directory to save data, default is data}
        --max_genome_length {maximum genome length to download, some of the genomes are huge and were not
                            used in jupyter notebook analysis}
'''


def url_to_json(url):
    r = requests.get(url)
    return r.json()

def make_page_url(page=None, page_size=None):
    if page is None and page_size is None:
        return "https://phagesdb.org/api/phages/"
    elif page is None:
        return "https://phagesdb.org/api/phages/?page_size={}".format(page_size)
    elif page_size is None:
        return "https://phagesdb.org/api/phages/?page={}".format(page)
    return "https://phagesdb.org/api/phages/?page={}&page_size={}".format(page, page_size)

def download_page(page=None, page_size=None):
    url = make_page_url(page=page, page_size=page_size)
    return url_to_json(url)

def download_phage_info(phage_name):
    url = "https://phagesdb.org/api/phages/{}/".format(phage_name)
    return url_to_json(url)

def download_sequenced_phages(page, page_size):
    #keep page = 1 and make page_size huge so can download in one go (unless memory is an issue)
    url = "https://phagesdb.org/api/sequenced_phages/?page={}&page_size={}".format(page, page_size)
    return url_to_json(url)

def download_url(url, pth):
    urllib.request.urlretrieve(url, os.path.join(pth, url.split("/")[-1]))

def main(pth, max_genome_length):
    if not os.path.isdir(pth):
        os.mkdir(pth)

    phage_pth = os.path.join(pth, "phages/")

    if not os.path.isdir(phage_pth):
        os.mkdir(phage_pth)

    print("Phage data will be saved at {}...".format(phage_pth))

    # Set an arbitrary upper limit that is bigger than the data
    phages_download = download_sequenced_phages(1, 10000)
    with open(os.path.join(pth, 'phage_data.txt'), 'w') as outfile:
        json.dump(phages_download, outfile)

    phages = phages_download['results']
    print("\tSucessfully downloaded info on {} phages".format(phages_download['count']))

    print("Downloading Phage Fasta Files...")

    for ii, pg in enumerate(phages):
        if ii % 500 == 0 and ii != 0:
            print("\t\t{}/{}".format(ii, len(phages)))

        if max_genome_length is None or int(pg['genome_length']) < max_genome_length:
            try:
                download_url(pg['fasta_file'], phage_pth)
            except:
                print("\t\t\tException on {}".format(pg['fasta_file']))

    print("\tSucesfully downloaded {} fasta files".format(len(os.listdir(phage_pth))))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--save_dir', default = "data/", help = "directory to store all the files")
    parser.add_argument('--max_genome_length', type = int, default = None,
                        help = "skips all phages w/ genomes longer than max_genome_length")
    args = parser.parse_args()
    main(args.save_dir, args.max_genome_length)
