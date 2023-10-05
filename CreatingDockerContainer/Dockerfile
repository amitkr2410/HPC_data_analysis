# Build from the official docker python base image, based on Debian

FROM python:3.6-stretch

# Install pre-reqs (commented ones are already in base image)

#RUN apt-get update && apt-get install -y emacs

#Install additional useful python packages

RUN pip install --upgrade pip  && pip install matplotlib numpy pandas scikit-learn

RUN cd /opt \
&& wget https://github.com/Kitware/CMake/releases/download/v3.17.0/cmake-3.17.0-Linux-x86_64.sh  \
&& echo "y" | bash ./cmake-3.17.0-Linux-x86_64.sh \
&& ln -s /opt/cmake-3.17.0-Linux-x86_64/bin/* /usr/local/bin

# Install Pythia8 
# ARG pythiaV="8307"
# RUN curl -SLk http://pythia.org/download/pythia83/pythia${pythiaV}.tgz \
# | tar -xvzC /usr/local \
# && cd /usr/local/pythia${pythiaV} \
# && ./configure --enable-shared --prefix=/usr/local  \
# && make -j2 \
# && make -j2 install


# Set environmental variables for cmake to know where things are (needed for SMASH, heppy)
ARG username=amit-docker
# ENV PYTHIA8DIR /usr/local/
# ENV PYTHIA8_ROOT_DIR /usr/local/
# ENV PATH $PATH:$PYTHIA8DIR/bin

# Set up a user group
ARG id=1234
RUN groupadd -g ${id} ${username} \
&& useradd --create-home --home-dir /home/${username} -u ${id} -g ${username} ${username}
USER ${username}
ENV HOME /home/${username}
WORKDIR ${HOME}

ENTRYPOINT /bin/bash
