# ChemTools Docker Configuration
# Based on miniconda2 for Python 2.7 + conda support
# Addresses GitHub Issue #42: Docker Update

FROM continuumio/miniconda2:latest

LABEL maintainer="ChemTools Dev Team <horton.chemtools@gmail.com>"
LABEL description="ChemTools - Interpretive Chemical Tools for Quantum Chemistry"
LABEL version="2.0"

# Set environment variables
ENV PYTHONDONTWRITEBYTECODE=1
ENV PYTHONUNBUFFERED=1

# Install system dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    git \
    git-lfs \
    && rm -rf /var/lib/apt/lists/* \
    && git lfs install

# Configure conda channels
RUN conda config --set always_yes yes && \
    conda config --add channels theochem && \
    conda config --add channels conda-forge

# Create conda environment and install dependencies
RUN conda create -n chemtools python=2.7 && \
    /bin/bash -c "source activate chemtools && \
    conda install -c theochem horton=2.1.0 && \
    conda clean -afy"

# Set up the conda environment to activate by default
ENV PATH="/opt/conda/envs/chemtools/bin:$PATH"
ENV CONDA_DEFAULT_ENV=chemtools

# Set working directory
WORKDIR /chemtools

# Copy project files
COPY . .

# Install ChemTools
RUN /bin/bash -c "source activate chemtools && pip install -e ."

# Create a non-root user for security
RUN useradd -m -s /bin/bash chemtools_user && \
    chown -R chemtools_user:chemtools_user /chemtools

USER chemtools_user

# Default command
ENTRYPOINT ["/bin/bash", "-c", "source activate chemtools && exec \"$@\"", "--"]
CMD ["chemtools", "--help"]
