#!/usr/bin/env python 
import ftplib
from ftplib import FTP
from io import StringIO
import mzIdentML_processor
ftp_server = ftplib.FTP('ftp.pride.ebi.ac.uk')
import pandas as pd
ftp_server.login()
path = '/pride/data/archive/2021/09/PXD027476/'
ftp_server.cwd(path)
ftp_server.dir()
filename = "README.txt"
with open(filename, "wb") as file:
    # Command for Downloading the file "RETR filename"
    ftp_server.retrbinary(f"RETR {filename}", file.write)
file = open(filename, "r")
for line in file:
    line_components = line.strip().split()
    if line_components[-2] == 'SEARCH':
        file_address = line_components[-3]
        if file_address.endswith('mzIdentML'):
            temp = file_address
            mz_name = temp.split('/')[9]
            file_path = str(path + mz_name + '/')
            with open(mz_name, "wb") as file:
                # Command for Downloading the file "RETR filename"
                ftp_server.retrbinary(f"RETR {mz_name}", file.write)
                # f = open(mz_name, "r")
                # mzIdentML_processor.processor(mz_name)
            mzIdentML_processor.processor(mz_name)
            break
