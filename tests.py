import unittest
import requests
import responses
from assessment import calculate_weight_of_each_crisp, find_crisprs_sequence


class TestMovieManagement(unittest.TestCase):

    test_off_target_summary_data = '{0: 1, 1: 2, 2: 3, 3: 4, 4: 5}'
    test_json_valid_response = {"species":"4","id":0,"off_target_summary":"{0: 1, 1: 2, 2: 3, 3: 4, 4: 5}","off_targets":[123456788,123456789]}
    test_json_invalid_response = {"error":"Error: 'ABCDEFGHTIKLMNOP' is an invalid argument for seq"}
    dna_sequence = "CTCAGCTAAGCACTCAGCTAAGGCAGCATCGATCTATCTTCATCTATTACTAGCGACTAGCATTATCATCG"

    def test_calculate_weight_of_each_crisp(self):
        calculated_weight = calculate_weight_of_each_crisp(self.test_off_target_summary_data)

        self.assertEqual(calculated_weight, 15)

    def test_find_crisprs_sequence_returns_valid_sequence(self):

        crispr_sequence = find_crisprs_sequence(self.dna_sequence)

        self.assertEqual(crispr_sequence, ['CTCAGCTAAGCACTCAGCTA'])

    @responses.activate
    def test_wge_stemcell_get_request_return_valid_response(self):

        responses.add(
            responses.GET,
            "https://wge.stemcell.sanger.ac.uk/api/off_targets_by_seq?seq=GATCGAGATAGTAATATGAT&species=Grch38&pam_right=false",
            json={"species":"4","id":0,"off_target_summary":"{0: 1, 1: 2, 2: 3, 3: 4, 4: 5}","off_targets":[123456788,123456789]},
            status=200,
        )

        mock_response = requests.get("https://wge.stemcell.sanger.ac.uk/api/off_targets_by_seq?seq=GATCGAGATAGTAATATGAT&species=Grch38&pam_right=false")

        self.assertDictEqual(mock_response.json(), self.test_json_valid_response)

    
    @responses.activate
    def test_wge_stemcell_get_request_return_invalid_response(self):

        responses.add(
            responses.GET,
            "https://wge.stemcell.sanger.ac.uk/api/off_targets_by_seq?seq=GATCGAGATAGTAATATGAT&species=Grch38&pam_right=false",
            json={"error":"Error: 'ABCDEFGHTIKLMNOP' is an invalid argument for seq"},
            status=400,
        )

        mock_response = requests.get("https://wge.stemcell.sanger.ac.uk/api/off_targets_by_seq?seq=GATCGAGATAGTAATATGAT&species=Grch38&pam_right=false")

        self.assertEqual(mock_response.status_code, 400)
        self.assertDictEqual(mock_response.json(), self.test_json_invalid_response)



if __name__ == '__main__':
    unittest.main()
