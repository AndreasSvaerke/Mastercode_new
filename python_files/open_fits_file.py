from pathlib import Path
from astropy.io import fits

def openfits(filepath):
    '''
    Takes filepath to folder with fits data and returns, 
    the names of the fits files, the data, and header in
    the same order as in the folder
    Returns: names, DATA, header
    '''
    data_file = filepath
    directory = Path(data_file)
    file_namesdata = [f.name for f in directory.iterdir() if f.is_file()]

    actudata = []
    for file_name in file_namesdata:
        actudata.append(file_name)
    
    DATA_data = []
    DATA_header0 = []
    #DATA_header = []

    for i in range(len(actudata)):
        with fits.open(data_file+'/'+actudata[i]) as hduldata:
            DATA_data.append(hduldata[0].data)
            DATA_header0.append(hduldata[0].header)
            #DATA_header.append(hduldata[1].header)

    return actudata, DATA_data, DATA_header0 #,DATA_header,

    