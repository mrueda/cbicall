import gdown

# Dictionary mapping the desired output filenames to their Google Drive file IDs
files = {
    'data.tar.gz.md5': '1jPi1YkQAxGaayKh_8XPmKPaq6HBAr4g6',
    'data.tar.gz.part-00': '14WR4RN3ohppYSJ1kpl5H1PCMivyEDGj4',
    'data.tar.gz.part-01': '1V1woKtshzi4w4tCEy-R2rvGAizi3woti',
    'data.tar.gz.part-02': '1DJOQBu3PqAk4nT6SAuE4PWQE8NmkMYHQ',
    'data.tar.gz.part-03': '1sXOKZiv4pZECQRAYiO4k2gjItqcDfrco'
}

# Loop over each file, build the download URL, and download it
for filename, file_id in files.items():
    url = f'https://drive.google.com/uc?export=download&id={file_id}'
    print(f"Downloading {filename}...")
    gdown.download(url, filename, quiet=False)
