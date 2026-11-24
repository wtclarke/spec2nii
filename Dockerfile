# syntax=docker/dockerfile:1

# Multi-stage build for spec2nii.
#   builder: installs git (for setuptools-scm versioning) and builds the
#            package and its dependencies into an isolated virtualenv.
#   runtime: slim Python image containing only the virtualenv, running as a
#            non-root user.

ARG PYTHON_VERSION=3.14

# ---------------------------------------------------------------------------
# Stage 1: builder
# ---------------------------------------------------------------------------
FROM python:${PYTHON_VERSION}-slim AS builder

ENV PIP_NO_CACHE_DIR=1 \
    PIP_DISABLE_PIP_VERSION_CHECK=1

# git is needed by setuptools-scm to derive the package version from tags.
RUN apt-get update \
 && apt-get install -y --no-install-recommends git \
 && rm -rf /var/lib/apt/lists/* \
 && git config --global --add safe.directory '*'

RUN python -m venv /opt/venv
ENV PATH="/opt/venv/bin:${PATH}"

WORKDIR /src
COPY . .

# Non-editable install into the venv; the venv is then copied to the runtime image.
RUN pip install --upgrade pip \
 && pip install .

# ---------------------------------------------------------------------------
# Stage 2: runtime
# ---------------------------------------------------------------------------
FROM python:${PYTHON_VERSION}-slim AS runtime

ARG UID=10001
ARG GID=10001

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PATH="/opt/venv/bin:${PATH}"

RUN groupadd --system --gid "${GID}" spec2nii \
 && useradd --system --uid "${UID}" --gid spec2nii --home-dir /home/spec2nii --create-home spec2nii

COPY --from=builder /opt/venv /opt/venv

# Working directory for mounted input/output data.
WORKDIR /data
RUN chown spec2nii:spec2nii /data

USER spec2nii

ENTRYPOINT ["spec2nii"]
CMD ["--help"]
