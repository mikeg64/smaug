# Use the latest ubuntu image as base for the new image
# ubuntu is the image name and latest is a tag that 
# references a particular version of the image.
# In this case latest is always the latest LTS image
# at the time of writing this, it's 16.04.
#FROM ubuntu:latest
FROM nvidia/cuda:10.0-devel-ubuntu16.04

#
RUN apt-get -y update
RUN apt-get -y upgrade
RUN apt-get install -y build-essential
RUN apt-get install -y gfortran
RUN apt-get install -y vim
RUN apt-get install -y nano
RUN apt-get install -y wget
RUN apt-get install -y curl

#RUN apt-get install -y firefox
# Run a system update to get it up to speed
# Then install python3 and pip
RUN apt-get update && apt-get install -y python2.7  \
    python-pip && apt-get install -y git-core
#RUN apt-get install -y python-tk

# Install jupyter
RUN pip2 install --upgrade pip
RUN pip2 install jupyter
RUN pip2 install matplotlib
RUN pip2 install pandas
RUN pip2 install -U scikit-learn
RUN pip2 install astropy
RUN pip2 install astropy-helpers
RUN pip2 install astroquery
RUN pip2 install aplpy
RUN pip2 install spectral-cube
RUN pip2 install reproject

RUN pip2 install sunpy[all]




# Add Tini. Tini operates as a process subreaper for jupyter. This prevents
# kernel crashes.
ENV TINI_VERSION v0.6.0
ADD https://github.com/krallin/tini/releases/download/${TINI_VERSION}/tini /usr/bin/tini
RUN chmod +x /usr/bin/tini
ENTRYPOINT ["/usr/bin/tini", "--"]




# Create a new system user
RUN useradd -ms /bin/bash jupyter

# Change to this new user
USER jupyter

# Set the container working directory to the user home folder
WORKDIR /home/jupyter
ENV HOME /home/jupyter

#RUN git clone https://github.com/littlelogking/virtual-astro
RUN git clone https://github.com/mikeg64/solar
RUN git clone https://github.com/mikeg64/smaug
RUN git clone https://github.com/mikeg64/sac_ccpforge
RUN git clone https://github.com/mikeg64/hermes

#pysac install failed
#RUN git clone https://github.com/SWAT-Sheffield/pysac
#RUN python pysac/setup.py -install

RUN ["mkdir", "notebooks"]
RUN ["mkdir", "notebooks/gallery_python"]
RUN ["mkdir", "notebooks/gallery_jupyter"]
COPY gallery_python/*.py notebooks/gallery_python/ 
COPY gallery_jupyter/*.ipynb notebooks/gallery_jupyter/

# Replace 1000 with your user / group id
#RUN export uid=1000 gid=1000 && \
#    echo "jupyter:x:${uid}:${gid}:Developer,,,:/home/jupyter:/bin/bash" >> /etc/passwd && \
#    echo "jupyter:x:${uid}:" >> /etc/group && \
#    echo "jupyter ALL=(ALL) NOPASSWD: ALL" > /etc/sudoers.d/jupyter && \
#    chmod 0440 /etc/sudoers.d/jupyter && \
#    chown ${uid}:${gid} -R /home/jupyter




#Jupyter
EXPOSE 8888
CMD ["jupyter", "notebook", "--port=8888", "--no-browser", "--ip=0.0.0.0"]


#to run this image
#running jupyter
#sudo docker run --rm -p 8888:8888 jupyter-sunpy



#docker run -ti --rm \
#       -e DISPLAY=$DISPLAY \
#       -v /tmp/.X11-unix:/tmp/.X11-unix \
#       firefox
#sudo docker run -it --rm -e DISPLAY=$DISPLAY -v /tmp/.X11-unix:/tmp/.X11-unix  jupyter-sunpy /bin/bash




