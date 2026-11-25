#!/usr/bin/env python3
"""Download a release archive, or use a locally downloaded copy, and extract it."""
import argparse
import hashlib
import json
from pathlib import Path
from urllib.request import urlretrieve
import zipfile


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('asset',choices=['models','predictions','force_constants','data'])
    ap.add_argument('--archive',type=Path,help='Use a local ZIP instead of downloading it')
    ap.add_argument('--dest',type=Path,default=Path('.'))
    args=ap.parse_args()
    info=json.loads((Path(__file__).resolve().parents[1]/'configs/assets.json').read_text())[args.asset]
    archive=args.archive
    if archive is None:
        if not info['url']:ap.error('No download URL configured. Supply --archive /path/to/'+info['filename'])
        archive=Path('downloads')/info['filename'];archive.parent.mkdir(parents=True,exist_ok=True)
        if not archive.exists():urlretrieve(info['url'],archive)
    digest=hashlib.sha256()
    with archive.open('rb') as f:
        for chunk in iter(lambda:f.read(2**20),b''):digest.update(chunk)
    if digest.hexdigest()!=info['sha256']:raise ValueError('Archive SHA256 does not match the release manifest.')
    args.dest.mkdir(parents=True,exist_ok=True);destination=args.dest.resolve()
    with zipfile.ZipFile(archive) as z:
        for name in z.namelist():
            if not (destination/name).resolve().is_relative_to(destination):
                raise ValueError('Archive path escapes the destination.')
        z.extractall(destination)
    print('Extracted',args.asset,'to',destination)


if __name__=='__main__':main()
