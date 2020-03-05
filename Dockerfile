FROM continuumio/anaconda3
COPY environment.yml .
COPY msiproc.py .
COPY peakselection.py .
COPY mir_parser.py .
RUN conda env create -f environment.yml
ENV PATH /opt/conda/envs/msiparse/bin:$PATH
RUN /bin/bash -c "source activate msiparse"
ENTRYPOINT ["python", "msiproc.py"]