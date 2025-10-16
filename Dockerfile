FROM python:3.14-slim

# Metadata
LABEL maintainer="Adriano Rutz <adafede@gmail.com>" \
      org.opencontainers.image.source="https://github.com/Adafede/wd-labels-to-iupac"

# Create non-root user with a real home directory
RUN useradd --no-log-init -r -u 1001 -g nogroup -d /home/nonroot -m nonroot

# Create application and cache directories
RUN mkdir -p /app /home/nonroot/.cache/uv && \
    chown -R nonroot:nogroup /app /home/nonroot

# Set environment variables
ENV UV_HOME=/app/.uv \
    HOME=/home/nonroot \
    PYTHONUNBUFFERED=1 \
    PIP_NO_CACHE_DIR=1

# Set working directory
WORKDIR /app

# Install system dependencies safely
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        libxi6 \
        libxrender1 \
        libxtst6 \
        openjdk-21-jre-headless && \
    rm -rf /var/lib/apt/lists/*

# Install uv globally
RUN --mount=type=cache,target=/root/.cache/pip \
    pip install --no-cache-dir uv

# Copy dependency file with proper ownership
COPY --chown=nonroot:nonroot pyproject.toml ./

# Install dependencies using uv with build cache
RUN --mount=type=cache,target=/app/.uv/cache \
    uv sync --cache-dir /app/.uv/cache

# Copy application code
COPY --chown=nonroot:nonroot wd_labels_to_iupac.py ./

# Use non-root user for security
USER nonroot

# Default command
CMD ["uv", "run", "wd_labels_to_iupac.py"]
