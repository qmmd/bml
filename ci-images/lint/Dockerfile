FROM ubuntu:xenial
COPY emacs.list /etc/apt/sources.list.d
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys 873503A090750CDAEB0754D93FF0E01EEAAFC9CD
RUN apt-get update
RUN apt-get install --assume-yes --no-install-recommends apt-utils
RUN apt-get install --assume-yes --no-install-recommends bundler
RUN apt-get install --assume-yes --no-install-recommends emacs26
RUN apt-get install --assume-yes --no-install-recommends git-core
RUN apt-get install --assume-yes --no-install-recommends indent
RUN apt-get install --assume-yes --no-install-recommends openssh-client
RUN apt-get install --assume-yes --no-install-recommends python-pip
RUN apt-get install --assume-yes --no-install-recommends python-setuptools
RUN apt-get install --assume-yes --no-install-recommends python-wheel
RUN pip install proselint
