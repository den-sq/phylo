# Build Stage
FROM rust AS build-stage

ADD . /usr/src/phylo
WORKDIR /usr/src/myapp

RUN cargo build --release

# Final Stage
FROM scratch

ARG GIT_COMMIT
ARG VERSION
LABEL REPO="https://github.com/HalfBrickInASock/phylo"
LABEL GIT_COMMIT=$GIT_COMMIT
LABEL VERSION=$VERSION

WORKDIR /usr/local/bin

COPY --from=build-stage /usr/src/phylo/bin/phylo /opt/phylo/bin/
RUN chmod +x /usr/local/bin/phylo

CMD /usr/local/bin/phylo
