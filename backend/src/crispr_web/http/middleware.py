"""ASGI middleware and direct-loopback request policy."""

from __future__ import annotations

import ipaddress
import logging

from fastapi import Request, status
from fastapi.responses import JSONResponse


LOGGER = logging.getLogger("crispr_web.api")


def _is_direct_loopback_request(request: Request) -> bool:
    if any(
        header in request.headers
        for header in ("x-forwarded-for", "x-real-ip", "forwarded")
    ):
        return False
    if request.client is None:
        return False
    try:
        address = ipaddress.ip_address(request.client.host)
    except ValueError:
        return False
    if isinstance(address, ipaddress.IPv6Address) and address.ipv4_mapped is not None:
        address = address.ipv4_mapped
    return address.is_loopback


class SecurityHeadersMiddleware:
    def __init__(self, app, *, hsts: bool):
        self.app = app
        self.hsts = hsts

    async def __call__(self, scope, receive, send):
        if scope["type"] != "http":
            await self.app(scope, receive, send)
            return

        async def send_with_headers(message):
            if message["type"] == "http.response.start":
                headers = list(message.get("headers", []))
                headers.extend(
                    [
                        (b"x-content-type-options", b"nosniff"),
                        (b"x-frame-options", b"DENY"),
                        (b"referrer-policy", b"no-referrer"),
                        (
                            b"permissions-policy",
                            b"camera=(), microphone=(), geolocation=()",
                        ),
                        (
                            b"content-security-policy",
                            b"default-src 'none'; frame-ancestors 'none'; base-uri 'none'",
                        ),
                        (b"cache-control", b"no-store"),
                    ]
                )
                if self.hsts:
                    headers.append(
                        (
                            b"strict-transport-security",
                            b"max-age=31536000; includeSubDomains",
                        )
                    )
                message["headers"] = headers
            await send(message)

        await self.app(scope, receive, send_with_headers)


class RequestLimitMiddleware:
    """Reject unbounded/chunked public submissions before body parsing."""

    def __init__(self, app, *, max_bytes: int):
        self.app = app
        self.max_bytes = max_bytes

    async def __call__(self, scope, receive, send):
        if scope["type"] != "http" or not (
            scope.get("method") == "POST" and scope.get("path") == "/api/v1/jobs"
        ):
            await self.app(scope, receive, send)
            return
        header_items = tuple(scope.get("headers", ()))
        raw_lengths = [
            value for key, value in header_items if key.lower() == b"content-length"
        ]
        has_transfer_encoding = any(
            key.lower() == b"transfer-encoding" for key, _value in header_items
        )
        if has_transfer_encoding or len(raw_lengths) > 1:
            response = JSONResponse(
                status_code=status.HTTP_400_BAD_REQUEST,
                content={"detail": "request framing is ambiguous"},
            )
            await response(scope, receive, send)
            return
        if not raw_lengths:
            response = JSONResponse(
                status_code=status.HTTP_411_LENGTH_REQUIRED,
                content={"detail": "Content-Length is required"},
            )
            await response(scope, receive, send)
            return
        raw_length = raw_lengths[0]
        try:
            content_length = int(raw_length)
        except ValueError:
            content_length = -1
        if content_length < 1:
            response = JSONResponse(
                status_code=status.HTTP_400_BAD_REQUEST,
                content={"detail": "request body is empty or has an invalid length"},
            )
            await response(scope, receive, send)
            return
        if content_length > self.max_bytes:
            response = JSONResponse(
                status_code=status.HTTP_413_CONTENT_TOO_LARGE,
                content={"detail": "request exceeds the public upload limit"},
            )
            await response(scope, receive, send)
            return
        messages = []
        received_bytes = 0
        while True:
            message = await receive()
            if message.get("type") == "http.disconnect":
                # The peer is gone. Do not invoke request parsing or generate a
                # noisy application error for a partial submission.
                return
            messages.append(message)
            if message.get("type") == "http.request":
                received_bytes += len(message.get("body", b""))
                if received_bytes > self.max_bytes:
                    response = JSONResponse(
                        status_code=status.HTTP_413_CONTENT_TOO_LARGE,
                        content={"detail": "request exceeds the public upload limit"},
                    )
                    await response(scope, receive, send)
                    return
                if not message.get("more_body", False):
                    break
        if received_bytes != content_length:
            response = JSONResponse(
                status_code=status.HTTP_400_BAD_REQUEST,
                content={"detail": "request body length does not match Content-Length"},
            )
            await response(scope, receive, send)
            return

        position = 0

        async def replay_receive():
            nonlocal position
            if position < len(messages):
                message = messages[position]
                position += 1
                return message
            return {"type": "http.disconnect"}

        await self.app(scope, replay_receive, send)


class JsonErrorMiddleware:
    """Return bounded JSON for unhandled HTTP failures before CORS decoration."""

    def __init__(self, app):
        self.app = app

    async def __call__(self, scope, receive, send):
        if scope["type"] != "http":
            await self.app(scope, receive, send)
            return
        response_started = False

        async def tracked_send(message):
            nonlocal response_started
            if message["type"] == "http.response.start":
                response_started = True
            await send(message)

        try:
            await self.app(scope, receive, tracked_send)
        except Exception:
            LOGGER.exception("unhandled API request failure")
            if response_started:
                raise
            response = JSONResponse(
                status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
                content={"detail": "the backend could not complete this request"},
            )
            await response(scope, receive, send)

