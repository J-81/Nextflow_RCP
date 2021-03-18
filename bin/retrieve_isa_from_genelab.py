#! /usr/bin/env python
from urllib.request import urlopen, quote, urlretrieve
from json import loads
from re import search
import argparse

import requests

def _parse_args():
  """ Parse command line args.
  """
  parser = argparse.ArgumentParser()
  parser.add_argument('--accession', metavar='GLDS-001', required=True,
                      help='GLDS accesion number')
  parser.add_argument('--alternate_url', action="store_true", default=False,
                      help='Use alternate url, fetched by api script')

  args = parser.parse_args()
  return args

# Function to pull metadata zip from GeneLab
# Credit to Kirill Grigorev
GENELAB_ROOT = "https://genelab-data.ndc.nasa.gov"
GLDS_URL_PREFIX = GENELAB_ROOT + "/genelab/data/study/data/"
FILELISTINGS_URL_PREFIX = GENELAB_ROOT + "/genelab/data/study/filelistings/"
ISA_ZIP_REGEX = r'.*_metadata_.*[_-]ISA\.zip$'

def read_json(url):
    with urlopen(url) as response:
        return loads(response.read().decode())

def get_isa(accession: str):
    glds_json = read_json(GLDS_URL_PREFIX + accession)
    try:
        _id = glds_json[0]["_id"]
    except (AssertionError, TypeError, KeyError, IndexError):
        raise ValueError("Malformed JSON?")
    isa_entries = [
        entry for entry in read_json(FILELISTINGS_URL_PREFIX + _id)
        if search(ISA_ZIP_REGEX, entry["file_name"])
    ]
    if len(isa_entries) == 0:
        raise ValueError("Unexpected: no ISAs found")
    elif len(isa_entries) > 1:
        raise ValueError("Unexpected: multiple files match the ISA regex")
    else:
        entry = isa_entries[0]
        version = entry["version"]
        url = GENELAB_ROOT + entry["remote_url"] + "?version={}".format(version)
        alt_url = (
            GENELAB_ROOT + "/genelab/static/media/dataset/" +
            quote(entry["file_name"]) + "?version={}".format(version)
        )
        return entry["file_name"], version, url, alt_url

def download_isa(accession: str, alternate_url: bool = False):
    """ Downloads isa for given accession number.

    :param accession: GLDS accession number, e.g. GLDS-194
    :param alternate_url: if true, uses alternative url, both alternate and default url should fetch the same file
    """
    print(f"Accessing GeneLab API for ISA file. Accesion: {accession}")
    filename ,_, url, alt_url  = get_isa(accession)
    print(f"Successfully retrieved ISA file location from API.")
    use_url = url if not alternate_url else alt_url
    print(f"Downloading from {use_url}. Alternative URL used: {alternate_url}")
    r = requests.get(use_url)
    with open(filename, "wb") as f:
        f.write(r.content)
    print(f"Finished downloading ISA file: {filename}")

if __name__ == "__main__":
    args = _parse_args()
    download_isa(args.accession, args.alternate_url)
