#!/usr/bin/env python

import gzip as gz
import json
import sys
import zipfile
from argparse import ArgumentParser
from io import BytesIO

import pandas as pd
import pycurl
from Bio.SeqIO import parse

SEARCH_URL = "https://files.jgi.doe.gov/mycocosm_file_list/"
DOWNLOAD_URL = "https://files-download.jgi.doe.gov/download_files/"
# Preferred file locations/names, in order, when several candidate files match a type
LOCATION_PRIORITY = ['Filtered Models ("best")', "All models"]
NAME_PRIORITY = ["GeneCatalog", "primary_alleles", "secondary_alleles", "all_"]


def read_token(path):
    """
    Read a JGI Data Portal session token from file, adding the 'Bearer' prefix if missing
    """
    with open(path) as fh:
        token = fh.read().strip()
    if not token.lower().startswith("bearer "):
        token = f"Bearer {token}"
    return token


def curl_request(url, cookie=None, headers=None, postfields=None):
    """
    Perform an HTTP request with pycurl and return the (status_code, response_body) tuple
    """
    buf = BytesIO()
    c = pycurl.Curl()
    c.setopt(c.URL, url)
    c.setopt(c.WRITEDATA, buf)
    c.setopt(c.VERBOSE, False)
    if cookie is not None:
        c.setopt(c.COOKIEFILE, cookie)
    if headers is not None:
        c.setopt(c.HTTPHEADER, headers)
    if postfields is not None:
        c.setopt(c.POSTFIELDS, postfields)
    c.perform()
    status = c.getinfo(c.RESPONSE_CODE)
    c.close()
    return status, buf.getvalue()


def get_json(portal, cookie):
    """
    Search the JGI Mycocosm file list for a portal and return the parsed JSON response
    """
    url = f"{SEARCH_URL}?organism={portal}#&api_version=2&a=false&h=false&d=asc&p=1&x=10&t=simple"
    status, body = curl_request(url, cookie=cookie)
    if status != 200:
        raise RuntimeError(f"Failed to fetch file list for {portal}: HTTP {status}")
    return json.loads(body)


def parse_files(json_dict):
    """
    Flatten the JGI file list JSON response into a DataFrame of file metadata indexed by file id

    :param json_dict: dict, parsed JSON response from get_json
    :return: (organism_id, DataFrame) tuple. organism_id identifies the organism when
             requesting downloads via download_file.
    """
    organism = json_dict["organisms"][0]
    records = {}
    for f in organism["files"]:
        meta = f["metadata"]
        records[f["_id"]] = {
            "name": f["file_name"],
            "type": f["file_type"],
            "md5": f["md5sum"],
            "format": meta.get("file_format"),
            "location": ";".join(meta["portal"]["display_location"]),
            "size": f["file_size"],
            "date": f["modified_date"],
        }
    file_df = pd.DataFrame(records).T
    file_df["dt"] = pd.to_datetime(file_df["date"])
    return organism["id"], file_df


def select_file(file_df, ft="transcripts"):
    """
    Select the file id and name of the best matching file of the given type ('transcripts' or
    'proteins') among a portal's files, preferring (in order): files under "Filtered Models
    (best)" over "All models", and within those, the canonical "GeneCatalog" file over
    primary/secondary allele or other variants.

    :param file_df: DataFrame, as returned by parse_files
    :param ft: str, 'transcripts' or 'proteins'
    :return: (file_id, file_name) tuple, or (None, None) if no matching file was found
    """
    target_type = ft.rstrip("s")
    is_type = file_df["type"].apply(lambda t: target_type in t)
    is_fasta = file_df["format"] == "fasta"
    is_gz = file_df["name"].str.endswith(".gz")
    not_deflines = ~file_df["name"].str.contains("deflines")
    candidates = file_df[is_type & is_fasta & is_gz & not_deflines]
    if candidates.empty:
        return None, None

    def rank(patterns, value):
        for i, pattern in enumerate(patterns):
            if pattern in value:
                return i
        return len(patterns)

    ranked = candidates.assign(
        _loc_rank=candidates["location"].apply(lambda v: rank(LOCATION_PRIORITY, v)),
        _name_rank=candidates["name"].apply(lambda v: rank(NAME_PRIORITY, v)),
        _size=candidates["size"].astype(int),
    ).sort_values(["_loc_rank", "_name_rank", "_size"])
    file_id = ranked.index[0]
    return file_id, candidates.loc[file_id, "name"]


def download_file(organism_id, file_id, filename, outfile, token):
    """
    Request and download a single file from JGI. The download API returns a zip archive
    containing a file manifest plus the requested file, so extract just the requested file.
    """
    payload = json.dumps({"ids": {organism_id: [file_id]}, "api_version": "2"})
    headers = [
        "accept: application/json",
        f"Authorization: {token}",
        "Content-Type: application/json",
    ]
    status, body = curl_request(DOWNLOAD_URL, headers=headers, postfields=payload)
    if status != 200:
        raise RuntimeError(f"Failed to download {filename}: HTTP {status}")
    with zipfile.ZipFile(BytesIO(body)) as zf:
        matches = [n for n in zf.namelist() if n.endswith(filename)]
        if not matches:
            raise RuntimeError(f"{filename} not found in downloaded archive")
        with zf.open(matches[0]) as fh_in, open(outfile, "wb") as fh_out:
            fh_out.write(fh_in.read())


def main(args):
    """
    Main function
    """
    token = read_token(args.token)
    json_dict = get_json(args.portal, args.cookie)
    organism_id, file_df = parse_files(json_dict)

    if args.outfile:
        file_id, filename = select_file(file_df, ft="transcripts")
        if file_id is None:
            sys.stderr.write(f"No transcripts file found for {args.portal}\n")
            sys.exit(1)
        sys.stderr.write(f"Downloading {filename} to {args.outfile}\n")
        download_file(organism_id, file_id, filename, args.outfile, token)

    if args.protein_out:
        file_id, filename = select_file(file_df, ft="proteins")
        if file_id is None:
            sys.stderr.write(f"No proteins file found for {args.portal}\n")
            sys.exit(1)
        sys.stderr.write(f"Downloading {filename} to {args.protein_out}\n")
        download_file(organism_id, file_id, filename, args.protein_out, token)
        if args.taxidmap:
            with (
                open(args.taxidmap, "w") as fhout,
                gz.open(args.protein_out, "rt") as fhin,
            ):
                fhout.writelines(
                    f"{record.id}\t{args.taxid}\n" for record in parse(fhin, "fasta")
                )


if __name__ == "__main__":
    parser = ArgumentParser(description="Download transcript files from JGI Mycocosm")
    parser.add_argument(
        "-p", "--portal", help="Portal shorthand name (e.g. Aaoar1)", required=True
    )
    parser.add_argument(
        "-c", "--cookie", help="Cookie file for JGI file search", required=True
    )
    parser.add_argument(
        "-t",
        "--token",
        help="File containing a JGI Data Portal session token, used for downloads "
        "(Avatar menu > Copy My Session Token on the JGI Data Portal)",
        required=True,
    )
    parser.add_argument(
        "-o",
        "--outfile",
        help="Output file name for transcript file. Omit this to skip download of transcript file",
    )
    parser.add_argument(
        "--protein_out",
        help="Attempt to find protein file for portal and write to file",
    )
    parser.add_argument(
        "--taxidmap", help="Output file with protein id to taxid mapping"
    )
    parser.add_argument("--taxid", help="Taxid of portal", type=int)
    args = parser.parse_args()
    if not args.outfile and not args.protein_out:
        sys.stderr.write("No output file specified\n")
        sys.exit(0)
    if args.taxidmap and not args.taxid:
        sys.stderr.write("No taxid specified for mapping\n")
        sys.exit(0)
    main(args)
