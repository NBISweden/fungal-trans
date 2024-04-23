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
    
def find_files(root):
    """
    Find files in the XML tree
    Searches for the latest transcript file in the folder "Filtered Models ("best")/Transcripts"
    """
    f = ""
    for child in root:
        if child.tag == "folder" and child.attrib["name"] == "Files":
            for item in child:
                if item.tag == "folder" and item.attrib["name"] == "Annotation":
                    for annot_item in item:
                        if annot_item.tag == "folder" and annot_item.attrib["name"] == 'Filtered Models ("best")':
                            for filtered_item in annot_item:
                                if filtered_item.tag == "folder" and filtered_item.attrib["name"] == "Transcripts":
                                    timestamp = False
                                    for transcript_item in filtered_item:
                                        if transcript_item.tag == "file":
                                            filename=transcript_item.attrib["filename"]
                                            if "transcripts" in filename and filename.endswith(".nt.fasta.gz"):
                                                url = transcript_item.attrib["url"]
                                                _timestamp = tz_convert(transcript_item.attrib["timestamp"])
                                                if timestamp and _timestamp>timestamp:
                                                    sys.stderr.write(f"Found newer file {filename} with timestamp {_timestamp}\n")
                                                    f = url
                                                elif not timestamp:
                                                    sys.stderr.write(f"Found file {filename} with timestamp {_timestamp}\n")
                                                    f = url
                                                    timestamp = _timestamp
                                    break
    return f

def main(args):
    """
    Main function
    """
    root = get_xml(args.portal, args.cookie)
    f = find_files(root)
    url = f"{args.base}{f}"
    if not args.outfile:
        outfile = os.path.basename(f)
    else:
        outfile = args.outfile
    sys.stderr.write(f"Downloading {os.path.basename(f)} to {outfile}\n")
    with open(outfile, 'wb') as fhout:
        c = pycurl.Curl()
        c.setopt(c.COOKIEFILE, args.cookie)
        c.setopt(c.URL, url)
        c.setopt(c.WRITEDATA, fhout)
        c.setopt(c.VERBOSE, False)
        c.perform()
        c.close()
    

if __name__ == '__main__':
    parser = ArgumentParser(description='Download all genomes from JGI Mycocosm')
    parser.add_argument('-p', '--portal', help='Portal shorthand name (e.g. Aaoar1)', required=True)
    parser.add_argument('-c', '--cookie', help="Cookie file for JGI", )
    parser.add_argument('-o', '--outfile', help='Output file name. If not given, the file will be saved in the current directory with the remote filename')
    parser.add_argument('-b', '--base', help='Base URL for JGI downloads (default: https://genome-downloads.jgi.doe.gov)', default="https://genome-downloads.jgi.doe.gov")
    args = parser.parse_args()
    main(args)