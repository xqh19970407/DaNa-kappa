#!/usr/bin/env python3
"""Fetch only the pinned upstream files needed for inference, into a local cache."""
import argparse
import base64
import json
from pathlib import Path
from urllib.request import Request,urlopen
from dana_kappa.upstream import UPSTREAM_COMMIT, FILES


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--dest',type=Path,default=Path.home()/'.cache/dana-kappa/aliegnn')
    args=ap.parse_args()
    prefix='https://api.github.com/repos/dannyzekunren/Dataset_thermoconductivity_pred/contents/'
    for name in FILES:
        dest=args.dest/name;dest.parent.mkdir(parents=True,exist_ok=True)
        if not dest.exists():
            request=Request(prefix+name+'?ref='+UPSTREAM_COMMIT,headers={'User-Agent':'DANA-kappa'})
            with urlopen(request,timeout=30) as response:
                record=json.load(response)
            dest.write_bytes(base64.b64decode(record['content']))
        print(dest)


if __name__=='__main__':main()
