FROM conda/miniconda3

LABEL authors="Matías Giménez <mgimenez@pasteur.edu.uy> & Ignacio Ferres <iferres@pasteur.edu.uy>"

# Copy conda recipe 
COPY plaSquid.yml /

# Install conda env
RUN conda env create -f /plaSquid.yml && \
    conda clean -afy

# Export  PATH
ENV PATH=/usr/local/envs/plaSquid/bin:$PATH
