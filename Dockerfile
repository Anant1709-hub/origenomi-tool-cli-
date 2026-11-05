FROM python:3.12.3


# Install NCBI BLAST
RUN apt-get update && \
    apt-get install -y ncbi-blast+ && \
    apt-get clean

# Install mummer (for nucmer)
RUN apt-get update && \
    apt-get install -y mummer && \
    apt-get clean

# Install progressiveMauve (system-level, unmaintained)
RUN apt-get update && \
    apt-get install -y progressivemauve && \
    apt-get clean



# Install pygenomeviz
RUN pip install --no-cache-dir pygenomeviz

RUN pip install --no-cache-dir origenomi

WORKDIR /work

# Expose origenomi CLI directly
ENTRYPOINT ["origenomi"]

# Default command 
CMD ["--help"]
