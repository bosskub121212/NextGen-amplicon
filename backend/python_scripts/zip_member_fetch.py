#!/usr/bin/env python3
"""Fetch single files out of a large remote ZIP without downloading the whole archive.

MIDORI2 ships each GenBank release as a few big ZIP bundles (every gene, ~0.4–0.5 GB);
we need one or two members. HTTP Range requests read the ZIP's central directory and
then just the bytes of the wanted member.

  python3 zip_member_fetch.py URL MEMBER OUT [MEMBER OUT ...]
  python3 zip_member_fetch.py --list URL

Standard library only. Exit 0 when every member was written, 1 otherwise.
"""
import io
import shutil
import sys
import urllib.request
import zipfile

UA = {"User-Agent": "Mozilla/5.0 (NextGen-Amplicon database downloader)"}


class HTTPRangeFile(io.RawIOBase):
    """Seekable read-only file over HTTP Range requests."""

    def __init__(self, url):
        self.url = url
        self.pos = 0
        req = urllib.request.Request(url, headers={**UA, "Range": "bytes=0-0"})
        with urllib.request.urlopen(req, timeout=60) as r:
            cr = r.headers.get("Content-Range")          # "bytes 0-0/12345"
            if not cr or "/" not in cr:
                raise OSError("server does not support range requests")
            self.size = int(cr.rsplit("/", 1)[1])
            self.url = r.geturl()                        # follow redirects once

    def readable(self):
        return True

    def seekable(self):
        return True

    def tell(self):
        return self.pos

    def seek(self, off, whence=0):
        self.pos = {0: off, 1: self.pos + off, 2: self.size + off}[whence]
        return self.pos

    def read(self, n=-1):
        if n is None or n < 0:
            n = self.size - self.pos
        if n == 0 or self.pos >= self.size:
            return b""
        end = min(self.pos + n, self.size) - 1
        req = urllib.request.Request(self.url, headers={**UA, "Range": f"bytes={self.pos}-{end}"})
        for attempt in range(4):
            try:
                with urllib.request.urlopen(req, timeout=120) as r:
                    data = r.read()
                break
            except OSError:
                if attempt == 3:
                    raise
        self.pos += len(data)
        return data

    def readinto(self, b):
        d = self.read(len(b))
        b[:len(d)] = d
        return len(d)


def main(argv):
    if len(argv) >= 2 and argv[0] == "--list":
        with zipfile.ZipFile(io.BufferedReader(HTTPRangeFile(argv[1]), 1 << 20)) as z:
            for i in z.infolist():
                print(f"{i.file_size:>12}  {i.filename}")
        return 0
    if len(argv) < 3 or len(argv[1:]) % 2:
        print(__doc__, file=sys.stderr)
        return 2
    url, pairs = argv[0], list(zip(argv[1::2], argv[2::2]))
    ok = True
    with zipfile.ZipFile(io.BufferedReader(HTTPRangeFile(url), 8 << 20)) as z:
        names = set(z.namelist())
        for member, out in pairs:
            if member not in names:
                print(f"  not in archive: {member}", file=sys.stderr)
                ok = False
                continue
            tmp = out + ".part"
            with z.open(member) as src, open(tmp, "wb") as dst:
                shutil.copyfileobj(src, dst, 8 << 20)
            shutil.move(tmp, out)
            print(f"  wrote {out}")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
