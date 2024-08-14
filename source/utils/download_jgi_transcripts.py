#!/usr/bin/env python

import xml.etree.ElementTree as ET
import pycurl
from datetime import datetime
import pytz
import tempfile
import os
from argparse import ArgumentParser
import sys


def tz_convert(timestamp):
    """
    Converts timestamp to UTC
    """
    timezones = {
        "PDT": "America/Los_Angeles",
        "PST": "America/Los_Angeles"
                 }
    tz = timestamp.split(" ")[-2]
    timestamp = " ".join(timestamp.split(" ")[0:-2]+timestamp.split(" ")[-1:])
    t = datetime.strptime(timestamp, "%a %b %d %H:%M:%S %Y")
    tz = pytz.timezone(timezones[tz])
    return tz.localize(t).astimezone(pytz.utc)


def get_xml(portal, cookie):
    """
    Get XML tree from JGI portal
    The xml file is stored in a temporary file and read into an ElementTree object
    """
    url=f"https://genome-downloads.jgi.doe.gov/portal/ext-api/downloads/get-directory?organism={portal}"
    with tempfile.TemporaryFile() as fp:
        c = pycurl.Curl()
        c.setopt(c.COOKIEFILE, cookie)
        c.setopt(c.URL, url)
        c.setopt(c.WRITEDATA, fp)
        c.setopt(c.VERBOSE, False)
        c.perform()
        c.close()
        fp.seek(0)
        tree = ET.fromstring(fp.read())
    return tree


def check_filtered(root):
    url = False
    timestamp = False
    files = []
    for _folder in root.findall("folder"):
        _foldername = _folder.attrib["name"]
        # exclude the Mycocosm folder
        if _foldername == "Mycocosm":
            continue
        _files=_folder.findall(".//file")
        for f in _files:
            _filename = f.attrib["filename"]
            if "filteredmodels" in _filename.lower() and "transcripts" in _filename and _filename.endswith(".gz"):
                files.append(f)
    for _file in files:
        _filename = _file.attrib["filename"]
        _url = _file.attrib["url"]
        _timestamp = tz_convert(_file.attrib["timestamp"])
        if not timestamp:
            sys.stderr.write(f"Found file {_filename} with timestamp {_timestamp}\n")
            url = _url
            timestamp = _timestamp
        elif timestamp and _timestamp>timestamp:
            sys.stderr.write(f"Found newer file {_filename} with timestamp {_timestamp}\n")
            url = _url
            timestamp = _timestamp
    return url


def check_filtered_folder(root):
    url = False
    timestamp = False
    files = []
    for folder in root.findall("folder"):
        foldername = folder.attrib["name"]
        if foldername == "Mycocosm":
            continue
        for subfolder in folder.findall(".//folder"):
            subfoldername = subfolder.attrib["name"]
            if 'Filtered Models ("best")' in subfoldername:
                sys.stderr.write(f"Found folder {foldername}/{subfoldername}\n")
                for _subfolder in subfolder.findall(".//folder"):
                    _subfoldername = _subfolder.attrib["name"]
                    if _subfoldername == "Transcripts":
                        sys.stderr.write(f"Found folder {foldername}/{subfoldername}/{_subfoldername}\n")
                        files = _subfolder.findall(".//file")
                        break
                if len(files) > 0:
                    break
    for file_element in files:
        _filename = file_element.attrib["filename"]
        if "transcripts" in _filename or "GeneModels" in _filename and _filename.endswith(".gz"):
            _url = file_element.attrib["url"]
            _timestamp = tz_convert(file_element.attrib["timestamp"])
            # if timestamp is not set, set it to the first file found
            if not timestamp:
                sys.stderr.write(f"Found file {_filename} with timestamp {_timestamp}\n")
                url = _url
                timestamp = _timestamp
            # if timestamp is set, check if the current file is newer
            elif timestamp and _timestamp>timestamp:
                sys.stderr.write(f"Found newer file {_filename} with timestamp {_timestamp}\n")
                url = _url
                timestamp = _timestamp
    return url


def find_mycocosm_file(root):
    md5 = False
    url = False
    for folder in root.findall("folder"):
        foldername = folder.attrib["name"]
        if foldername == "Mycocosm":
            for subfolder in folder.findall(".//folder"):
                subfoldername = subfolder.attrib["name"]
                if 'Filtered Models ("best")' in subfoldername:
                    sys.stderr.write(f"Found folder {foldername}/{subfoldername}\n")
                    for _subfolder in subfolder.findall(".//folder"):
                        _subfoldername = _subfolder.attrib["name"]
                        if _subfoldername == "Transcripts":
                            sys.stderr.write(f"Found folder {foldername}/{subfoldername}/{_subfoldername}\n")
                            for f in _subfolder.findall(".//file"):
                                filename = f.attrib["filename"]
                                if "transcripts" in filename or "GeneModels" in filename or filename.endswith("fasta.gz") and filename.endswith(".gz"):
                                    md5 = f.attrib["md5"]
                                    break
                    if md5:
                        break
    if not md5:
        return False
    files = []
    for folder in root.findall("folder"):
        foldername = folder.attrib["name"]
        if foldername == "Mycocosm":
            continue
        files+=folder.findall(".//file")
    for f in files:
        if "md5" in f.attrib.keys() and f.attrib["md5"] == md5:
            url = f.attrib["url"]
            timestamp = tz_convert(f.attrib["timestamp"])
            sys.stderr.write(f"Found file by md5 checksum {f.attrib['filename']} with timestamp {timestamp}\n")
    return url

def find_smallest(root):
    url = False
    size = False
    sizeinbytes = False
    for folder in root.findall("folder"):
        foldername = folder.attrib["name"]
        if foldername == "Mycocosm":
            continue
        for f in folder.findall(".//file"):
            filename = f.attrib["filename"]
            _size = f.attrib["size"]
            _sizeinbytes = f.attrib["sizeInBytes"]
            if "transcripts" in filename and filename.endswith(".gz"):
                if not sizeinbytes:
                    sizeinbytes = int(_sizeinbytes)
                    size = _size
                    url = f.attrib["url"]
                    sys.stderr.write(f"Found file {filename} with size {size}\n")
                elif sizeinbytes and int(_sizeinbytes)<sizeinbytes:
                    sizeinbytes = int(_sizeinbytes)
                    url = f.attrib["url"]
                    sizeinbytes = _sizeinbytes
                    sys.stderr.write(f"Found smaller file {filename} with size {size}\n")
    return url

def find_files(root):
    """
    Find files in the XML tree
    Searches for the latest transcript file in the folder "Filtered Models ("best")/Transcripts"
    """
    url = False
    # check for files matching 'FilteredModels3.transcripts.fasta.gz'
    sys.stderr.write("Searching for FilteredModels transcript file\n")
    url = check_filtered(root)
    if url:
        return url
    sys.stderr.write("No FilteredModels file found\n")
    # if no url found, check folders to locate the Filtered Models ("best")/Transcripts folder
    sys.stderr.write('Searching for Filtered Models ("best")/Transcripts folder \n')
    url = check_filtered_folder(root)
    if url:
        return url
    # if no such folder exists with file, try searching for a filtered transcripts file under Mycocosm subfolder
    sys.stderr.write('Filtered Models ("best")/Transcripts folder found\n')
    sys.stderr.write("Searching for Mycocosm transcript file\n")
    url = find_mycocosm_file(root)
    if url:
        return url
    # as a last resort, download the smallest transcript file found
    sys.stderr.write("No file found, falling back to finding smallest transcript file\n")
    url = find_smallest(root)
    return url


def write_file(url, outfile, cookie):
    with open(outfile, 'wb') as fhout:
        c = pycurl.Curl()
        c.setopt(c.COOKIEFILE, cookie)
        c.setopt(c.URL, url)
        c.setopt(c.WRITEDATA, fhout)
        c.setopt(c.VERBOSE, True)
        c.perform()
        c.close()


def main(args):
    """
    Main function
    """
    root = get_xml(args.portal, args.cookie)
    url = find_files(root)
    if not url:
        sys.stderr.write("No files found\n")
        url=""
    url = f"{args.base}{url}"
    if not args.outfile:
        outfile = os.path.basename(url)
    else:
        outfile = args.outfile
    sys.stderr.write(f"Downloading {url} to {outfile}\n")
    write_file(url, outfile, args.cookie)
    

if __name__ == '__main__':
    parser = ArgumentParser(description='Download all genomes from JGI Mycocosm')
    parser.add_argument('-p', '--portal', help='Portal shorthand name (e.g. Aaoar1)', required=True)
    parser.add_argument('-c', '--cookie', help="Cookie file for JGI", )
    parser.add_argument('-o', '--outfile', help='Output file name. If not given, the file will be saved in the current directory with the remote filename')
    parser.add_argument('-b', '--base', help='Base URL for JGI downloads (default: https://genome-downloads.jgi.doe.gov)', default="https://genome-downloads.jgi.doe.gov")
    args = parser.parse_args()
    main(args)