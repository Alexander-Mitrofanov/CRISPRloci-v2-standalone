"""Exercise the actual Nginx maps: preflights must not spend download tokens."""

from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
import re
import shutil
import socket
import subprocess
import threading
import time
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen

import pytest


NGINX = shutil.which("nginx")
pytestmark = pytest.mark.skipif(NGINX is None, reason="Nginx executable unavailable")
ROOT = Path(__file__).resolve().parents[1]


class Upstream(BaseHTTPRequestHandler):
    def do_GET(self):
        self.send_response(200)
        self.end_headers()

    do_HEAD = do_GET

    def do_OPTIONS(self):
        self.send_response(204)
        self.end_headers()

    def log_message(self, *_args):
        pass


def status(origin, route, method):
    try:
        with urlopen(Request(origin + route, method=method), timeout=2) as response:
            return response.status
    except HTTPError as error:
        return error.code


@pytest.mark.parametrize("method", ["GET", "HEAD"])
def test_preflights_do_not_consume_but_downloads_still_exhaust_budget(tmp_path, method):
    template = (ROOT / "deploy/nginx/crispr-web-api.conf.template").read_text()
    maps = re.findall(
        r"^map [^\n]+ \$crispr_download_(?:addr|global) \{.*?^\}",
        template, re.MULTILINE | re.DOTALL,
    )
    assert len(maps) == 2
    # Only the client address transport differs from the deployed PROXY listener.
    maps = "\n".join(maps).replace("$proxy_protocol_addr", "$remote_addr")
    with socket.socket() as listener:
        listener.bind(("127.0.0.1", 0))
        port = listener.getsockname()[1]
    upstream = ThreadingHTTPServer(("127.0.0.1", 0), Upstream)
    thread = threading.Thread(target=upstream.serve_forever, daemon=True)
    thread.start()
    config = tmp_path / "nginx.conf"
    config.write_text(f"""
pid {tmp_path}/nginx.pid;
error_log {tmp_path}/error.log error;
events {{ worker_connections 64; }}
http {{
    access_log off;
    client_body_temp_path {tmp_path}/body;
    proxy_temp_path {tmp_path}/proxy;
    fastcgi_temp_path {tmp_path}/fastcgi;
    uwsgi_temp_path {tmp_path}/uwsgi;
    scgi_temp_path {tmp_path}/scgi;
    {maps}
    limit_req_zone $crispr_download_addr zone=downloads:1m rate=1r/m;
    limit_req_zone $crispr_download_global zone=global_downloads:1m rate=1r/m;
    server {{
        listen 127.0.0.1:{port};
        server_name localhost;
        limit_req_status 429;
        limit_req zone=downloads burst=1 nodelay;
        limit_req zone=global_downloads burst=1 nodelay;
        location / {{ proxy_pass http://127.0.0.1:{upstream.server_port}; }}
    }}
}}
""")
    process = None
    try:
        process = subprocess.Popen(
            [NGINX, "-p", str(tmp_path), "-c", str(config),
             "-g", "daemon off; master_process off;"],
            stdout=subprocess.DEVNULL, stderr=subprocess.PIPE,
        )
        origin = f"http://127.0.0.1:{port}"
        deadline = time.monotonic() + 5
        while True:
            assert process.poll() is None, process.communicate()[1].decode()
            try:
                assert status(origin, "/", "GET") == 200
                break
            except (URLError, TimeoutError):
                assert time.monotonic() < deadline
                time.sleep(0.02)
        base = "/api/v1/jobs/" + "a" * 32
        routes = [base + "/result", base + "/artifacts/" + "b" * 24]
        for route in routes * 2:
            assert status(origin, route, "OPTIONS") == 204
        for route in routes:
            assert status(origin, route, method) == 200
        assert status(origin, routes[0], method) == 429
    finally:
        if process is not None:
            process.terminate()
            process.communicate(timeout=5)
        upstream.shutdown()
        upstream.server_close()
        thread.join(timeout=5)
