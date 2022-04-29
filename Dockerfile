# Use offical r image with r version algorithm was tested on
FROM r-base:3.6.2

WORKDIR /Bayesian_LSP

COPY . /Bayesian_LSP

# Install dependencies
RUN apt-get update && \
    apt install -y jags && \
    Rscript ./requirements/requirements.R

# Make directory for data just incase tests need outputs
CMD mkdir /data

# tests
# CMD [ "./test.sh" ]

