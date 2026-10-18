FROM python:3.12-slim-bookworm AS build
RUN apt-get update && apt-get install -y --no-install-recommends g++ make git && rm -rf /var/lib/apt/lists/*
WORKDIR /src
COPY . .
RUN pip install --no-cache-dir -r requirements.lock build && cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release -DHAVE_SSE4_1=1 -DHAVE_AVX2=0 && cmake --build build --parallel 4 && ctest --test-dir build --output-on-failure && python -m build
FROM python:3.12-slim-bookworm
RUN apt-get update && apt-get install -y --no-install-recommends libgomp1 && rm -rf /var/lib/apt/lists/*
COPY --from=build /src/build/src/spacepharer-native /usr/local/bin/spacepharer-native
COPY --from=build /src/dist/*.whl /tmp/wheels/
RUN pip install --no-cache-dir /tmp/wheels/*.whl && rm -rf /tmp/wheels
ENTRYPOINT ["spacepharer"]
