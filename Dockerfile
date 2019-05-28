from ubuntu:19.04

RUN apt-get update && apt-get install -y \
    build-essential libeigen3-dev
