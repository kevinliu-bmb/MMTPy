import os
import requests

class FetchDiet:
    def __init__(self):
        self.base_url = 'https://www.vmh.life/_api/dietflux/'
        self.base_path = 'data/diet'

    def fetch_unique_diets(self):
        url = self.base_url
        diets = set()
        while url:
            response = requests.get(url)
            if response.status_code == 200:
                data = response.json()
                for entry in data['results']:
                    diets.add(entry['diet'])
                url = data.get('next')
            else:
                print("Failed to fetch data: Status code", response.status_code)
                return []
        return list(diets)

    def download_diet_flux_csv(self, diet, file_name='diet_flux.csv'):
        full_path = os.path.join(self.base_path, f'{diet.replace(" ", "_")}_{file_name}')

        if not os.path.exists(self.base_path):
            os.makedirs(self.base_path)

        url = f'{self.base_url}?diet={diet}&format=pcsv'
        response = requests.get(url)
        if response.status_code == 200:
            with open(full_path, 'wb') as file:
                file.write(response.content)
            print(f"File '{full_path}' has been downloaded.")
        else:
            print("Failed to download file: Status code", response.status_code)

