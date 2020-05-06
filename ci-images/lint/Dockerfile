FROM ubuntu:xenial

COPY emacs.list /etc/apt/sources.list.d
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys 873503A090750CDAEB0754D93FF0E01EEAAFC9CD

RUN apt-get update
RUN ln -fs /usr/share/zoneinfo/UTC /etc/localtime
RUN apt-get install --assume-yes tzdata
RUN DEBIAN_FRONTEND=noninteractive dpkg-reconfigure --frontend noninteractive tzdata

RUN apt-get install --assume-yes --no-install-recommends \
        apt-utils \
        bundler \
        emacs26 \
        git-core \
        indent \
        openssh-client \
        python-pip \
        python-setuptools \
        python-wheel

RUN pip install proselint
RUN groupadd --system ci
RUN useradd --system --gid ci --create-home ci
WORKDIR /home/ci
USER ci
