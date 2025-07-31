from datetime import timedelta
import pandas as pd
import sys

if len(sys.argv) < 3:
    print("Usage: python generate_urls.py <ensemble_number> <output_folder>")
    sys.exit(1)

ne = int(sys.argv[1])
output_dir = sys.argv[2]
# 'https://data.nssl.noaa.gov/thredds/dodsC/WRDD/Winter-PHI/20200130/0000/ENS_MEM_01/wrfout_d01_2020-01-31_00:20:00'

date_str = '20200130'

#initializations = ['2000', '2100', '2200', '2300', '0000', '0100', '0200', '0300']


initializations = ['0100']
ensembles = [f'ENS_MEM_{ne:02d}']
base_url = 'https://data.nssl.noaa.gov/thredds/fileServer/WRDD/Winter-PHI'

file_start = '/wrfout_d01_'

# Build base path
base_path = f"{base_url}/{date_str}/"

# Output list of filenames
file_list = []
for hour in initializations:
    path = base_path + hour + '/'

    for ens in ensembles:
        next_path = path + ens
        first_date = pd.to_datetime(date_str+hour, format='%Y%m%d%H%M')
        if hour[:1] == '0':
            first_date = first_date + timedelta(days=1)

        for file in range(73):
            new_date_str = first_date.strftime('%Y-%m-%d_%H:%M')
            file_name = next_path + file_start + new_date_str + ':00'
            print(file_name)
            file_list.append(file_name)

            first_date = first_date + timedelta(minutes=5)

output_path = f'{output_dir}/wofs_file_list.txt'

with open(output_path, 'w') as f:
    for url in file_list:
        f.write(url + '\n')

print("✅ File list saved to '{output_path}'")
