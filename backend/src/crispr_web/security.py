"""Capability tokens and trustworthy client-address extraction."""

from __future__ import annotations

import hashlib
import hmac
import ipaddress
import secrets
from typing import Iterable

from fastapi import Request


def new_access_token() -> str:
    """Return exactly 256 random bits encoded as URL-safe text."""

    return secrets.token_urlsafe(32)


def token_digest(token: str) -> str:
    return hashlib.sha256(token.encode("ascii", errors="strict")).hexdigest()


def token_matches(token: str, expected_digest: str) -> bool:
    try:
        observed = token_digest(token)
    except (UnicodeError, AttributeError):
        return False
    return hmac.compare_digest(observed, expected_digest)


def _normalized_address(
    value: str,
) -> ipaddress.IPv4Address | ipaddress.IPv6Address:
    """Parse one unscoped address and collapse IPv4-mapped IPv6 to IPv4."""

    address = ipaddress.ip_address(value)
    if isinstance(address, ipaddress.IPv6Address):
        if address.scope_id is not None:
            raise ValueError("scoped IPv6 addresses are not valid client identities")
        if address.ipv4_mapped is not None:
            return address.ipv4_mapped
    return address


def client_rate_identity(client_address: str) -> str:
    """Return the network identity used for privacy-preserving admission limits.

    IPv4 clients retain a canonical host identity. IPv6 clients are grouped by
    their /64 allocation so temporary privacy addresses cannot evade the active
    job or submission limits. The queue persists only a keyed digest of this
    value, never the raw address or network.
    """

    address = _normalized_address(client_address)
    if isinstance(address, ipaddress.IPv4Address):
        return address.compressed
    network = ipaddress.IPv6Network((address, 64), strict=False)
    return network.with_prefixlen


def client_ip(
    request: Request,
    trusted_proxy_networks: Iterable[ipaddress.IPv4Network | ipaddress.IPv6Network],
) -> str:
    """Resolve an IP without trusting forwarding headers from public peers.

    Uvicorn must be started with ``--no-proxy-headers`` so the immediate peer
    remains available here.  When that peer is an explicitly configured local
    proxy, the chain is evaluated right-to-left and only explicitly trusted hops are removed.
    """

    peer_text = request.client.host if request.client is not None else "0.0.0.0"
    try:
        peer = _normalized_address(peer_text)
    except ValueError:
        peer = ipaddress.IPv4Address("0.0.0.0")
    forwarded = request.headers.get("x-forwarded-for", "")
    if not forwarded:
        return peer.compressed
    chain = []
    for value in forwarded.split(","):
        try:
            chain.append(_normalized_address(value.strip()))
        except ValueError:
            return peer.compressed
    current = peer
    for candidate in reversed(chain):
        if not any(current in network for network in trusted_proxy_networks):
            break
        current = candidate
    return current.compressed


def bearer_token(authorization: str | None) -> str | None:
    if not authorization:
        return None
    scheme, separator, token = authorization.partition(" ")
    if not separator or scheme.lower() != "bearer" or not token.strip():
        return None
    value = token.strip()
    if len(value) > 256 or any(character.isspace() for character in value):
        return None
    return value
