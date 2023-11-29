import unittest
from unittest.mock import Mock, patch, mock_open
import sys

sys.path.append('src/utils')

from fetch_diet import FetchDiet

class TestFetchDiet(unittest.TestCase):
    @patch('requests.get')
    def test_fetch_unique_diets(self, mock_get):
        fetch_diet = FetchDiet()
        mock_response = Mock()
        mock_response.status_code = 200
        mock_response.json.return_value = {'results': [{'diet': 'Diet A'}, {'diet': 'Diet B'}]}
        mock_get.return_value = mock_response
    
        diets = fetch_diet.fetch_unique_diets()
        self.assertCountEqual(diets, ['Diet A', 'Diet B'])

    @patch('os.path.exists', return_value=False)
    @patch('os.makedirs')
    @patch('requests.get')
    @patch('builtins.open', new_callable=mock_open)
    def test_download_diet_flux_csv(self, mock_open, mock_get, mock_makedirs, mock_exists):
        fetch_diet = FetchDiet()
        mock_response = Mock()
        mock_response.status_code = 200
        mock_response.content = b'Dummy CSV Content'
        mock_get.return_value = mock_response
    
        fetch_diet.download_diet_flux_csv('Diet A')
    
        mock_makedirs.assert_called_once_with('data/diet')
        mock_open.assert_called_once_with('data/diet/Diet_A_diet_flux.csv', 'wb')
        mock_open().write.assert_called_once_with(b'Dummy CSV Content')

if __name__ == "__main__":
    unittest.main()