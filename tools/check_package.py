"""Validate the distributable file set, local HTML links and source archives."""
from html.parser import HTMLParser
from pathlib import Path
from urllib.parse import unquote, urlsplit
import gzip
import hashlib
import json
import subprocess
import zipfile

ROOT = Path(__file__).resolve().parents[1]


class Links(HTMLParser):
    def __init__(self):
        super().__init__()
        self.links = []

    def handle_starttag(self, tag, attrs):
        self.links.extend(value for key, value in attrs if key in ('href', 'src') and value)


def main():
    names = subprocess.check_output(['git', 'ls-files', '--cached', '--others', '--exclude-standard'], cwd=ROOT, text=True).splitlines()
    names = sorted({name for name in names if (ROOT / name).is_file()})
    total = sum((ROOT / name).stat().st_size for name in names)
    assert total < 1_000_000_000, total
    assert all((ROOT / name).stat().st_size < 100 * 1024**2 for name in names)
    assert not any(name.startswith(('.local/', 'tmp/', '.claude/', 'publication/')) for name in names)
    targets = set(names)
    checked = 0
    for name in names:
        if not name.endswith('.html') or name.endswith('.template.html'):
            continue
        parser = Links()
        parser.feed((ROOT / name).read_text(encoding='utf8'))
        for link in parser.links:
            url = urlsplit(link)
            if url.scheme or url.netloc or not url.path:
                continue
            path = (ROOT / name).parent / unquote(url.path)
            target = path.resolve().relative_to(ROOT).as_posix()
            assert target in targets, (name, link)
            checked += 1
    for item in json.loads((ROOT / 'distribution.json').read_text())['compressed_sources']:
        raw = (ROOT / item['archive']).read_bytes()
        assert hashlib.sha256(raw).hexdigest() == item['archive_sha256']
        restored = gzip.decompress(raw)
        assert len(restored) == item['bytes'] and hashlib.sha256(restored).hexdigest() == item['sha256']
    directory = ROOT / 'research/symmetrizable'
    manifest = json.loads((directory / 'source-manifest.json').read_text())
    archive = directory / 'certificate-sources.zip'
    assert hashlib.sha256(archive.read_bytes()).hexdigest() == manifest['archive_sha256']
    with zipfile.ZipFile(archive) as sources:
        for item in manifest['files']:
            raw = sources.read(item['path'])
            assert raw == (ROOT / item['path']).read_bytes(), item['path']
            assert hashlib.sha256(raw).hexdigest() == item['sha256'], item['path']
        for item in manifest['historical_enumerators']:
            assert hashlib.sha256(sources.read(item['archive_member'])).hexdigest() == item['source_sha256']
    print(f'PASS: {len(names)} distributable files ({total / 1024**2:.2f} MiB), {checked} local HTML links, compressed inputs and archived sources.')


if __name__ == '__main__':
    main()
