import logging
import requests
from requests.exceptions import HTTPError, ConnectTimeout, ConnectionError
import csv


def find_crisprs_sequence(dna_sequence):

    """
        This function takes dna_sequence as an argument and returns the crisp sequence based on given search 
        criteria.
    """
    crisprs = []
    sequence_length = len(dna_sequence)
    guide_length = 20 
    pam_sequence_length = 3
    pam_sequence = ["AGG", "TGG"]

    for i in range(sequence_length - guide_length - pam_sequence_length + 1):
        crispr_candidate = dna_sequence[i:i + guide_length]
        pam_candidate = dna_sequence[i + guide_length:i + guide_length +
                                        pam_sequence_length]
        if pam_candidate in pam_sequence:
            crisprs.append(crispr_candidate)

    return crisprs


def calculate_weight_of_each_crisp(off_target_summary):
    """
        Calculate total weight of each crisprs from off target summary.
    """
    return sum(eval(off_target_summary).values())

    
def main():
    try:
        tsv_file = open("Sequences.tsv")
        data = csv.DictReader(tsv_file, delimiter="\t")
    except (FileNotFoundError, AttributeError) as exc:
        logging.error(f"Error occurred while reading TSV file with Exception {exc}")

    result = []

    if not data:
        logging.debug("No data found in current file.")
        return
    
    for dna_sequence in data:
        crisp_sequences = find_crisprs_sequence(dna_sequence['Seq'])
        if len(crisp_sequences): # Skip empty sequence
            result.append({"ID": dna_sequence['ID '], "crisp_seq": crisp_sequences})

    off_target_summary = {}

    for crisp_seq in result:
        for seq in crisp_seq.get('crisp_seq'):
            url = f"https://wge.stemcell.sanger.ac.uk/api/off_targets_by_seq?seq={seq}&species=Grch38&pam_right=true"
            try:
                response = requests.get(url)
                weight = calculate_weight_of_each_crisp(response.json()["off_target_summary"])
                off_target_summary.update({seq: weight})
            except (HTTPError, ConnectTimeout, ConnectionError) as exc:
                logging.error(f"An Error occurred while requesting to url {url} /n with exception {exc}")
            

    sorted_crisp_seq = sorted(off_target_summary.items(), key=lambda x: x[1], reverse=True)

    for i, seq in enumerate(sorted_crisp_seq, 1):
        print(f"- Crisper Sequence {i}: {seq[0]} and it's Weight: {seq[1]}")


  
if __name__ == "__main__":
    main()
  