running and building docker image with a gpu

To build this using docker needs nvidia version of docker installed as follows



Built using 
sudo docker build . -t sac-smaug-docker

Run docker using
sudo docker run --runtime=nvidia --rm -p 8888:8888 sac-smaug-docker


Run bash terminal using docker run -it --rm mikeg64/sac-smaug-docker /bin/bash

create a volume on the host on windows this may be in C:\ProgramData\Docker
docker volume create --name sunpy-vol --driver local

To access the volume using the created volume use
docker run -v C:\ProgramData\Docker\sunpy-vol:/home/jupyter/solar/docker --runtime=nvidia  --rm -p 8888:8888 mikeg6
4/sac-smaug-docker



https://marmelab.com/blog/2018/03/21/using-nvidia-gpu-within-docker-container.html
Using NVIDIA GPU within Docker Containers
Jonathan Petitcolas
Jonathan PetitcolasMarch 21, 2018
#devops#ai

Diving into machine learning requires some computation power, mainly brought by GPUs. But I'm reluctant to install new software stacks on my laptop - I prefer installing them in Docker containers, to avoid polluting other programs, and to be able to share the results with my coworkers. That means I have to configure Docker to use my GPU. This is the story of how I managed to do it, in about half a day.
Introduction to NVIDIA Docker

I'm used to using Docker for all my projects at marmelab. It allows to setup easily even the most complex infrastructures, without polluting the local system. However, as image processing generally requires a GPU for better performances, the first question is: can Docker handle GPUs?

Looking for an answer to this question leads me to the nvidia-docker repository, described in a concise and effective way as:

    Build and run Docker containers leveraging NVIDIA GPUs

Fortunately, I have an NVIDIA graphic card on my laptop. NVIDIA engineers found a way to share GPU drivers from host to containers, without having them installed on each container individually.

Enable NVIDIA GPU within Docker containers

GPUs on container would be the host container ones. Looks promising. Let's give it a try!
Installing CUDA on Host

CUDA is a parallel computing platform allowing to use GPU for general purpose processing. It is strongly recommended when dealing with machine learning, an important resource consuming task.
Does our Graphical Card supports CUDA?

The first step is to identify precisely the model of my graphical card. This is done easily on Linux using the lspci util:

lspci | grep VGA

00:02.0 VGA compatible controller: Intel Corporation [...] (rev 09)
01:00.0 VGA compatible controller: NVIDIA [...] [GeForce GT 650M] (rev a1)

So, I have a GT650M. Checking on NVIDIA commercial website, this card has CUDA support. Great so far!


Having Result = PASS shows that my CUDA installation is fully operational!
Installing NVIDIA Docker

Installing nvidia-docker is far easier than installing CUDA. First, I need to add the nvidia-docker dependencies:

curl -s -L https://nvidia.github.io/nvidia-docker/gpgkey | sudo apt-key add -
curl -s -L https://nvidia.github.io/nvidia-docker/ubuntu16.04/amd64/nvidia-docker.list | \
    sudo tee /etc/apt/sources.list.d/nvidia-docker.list

sudo apt-get update

Before installing nvidia-docker2 utility, I need to ensure that I use docker-ce, the latest official Docker release. Based on official documentation, here is the process to follow:

# remove all previous Docker versions
sudo apt-get remove docker docker-engine docker.io

# add Docker official GPG key
curl -fsSL https://download.docker.com/linux/ubuntu/gpg | sudo apt-key add -

# Add Docker repository (for Ubuntu Xenial)
sudo add-apt-repository \
    "deb [arch=amd64] https://download.docker.com/linux/ubuntu xenial stable"

sudo apt-get update
sudo apt install docker-ce

Now that I have the last version of Docker, I can install nvidia-docker:

# Install nvidia-docker2 and reload the Docker daemon configuration
sudo apt-get install -y nvidia-docker2
sudo pkill -SIGHUP dockerd

I now have access to a Docker nvidia runtime, which embeds my GPU in a container. I can use it with any Docker container.

Let's ensure everything work as expected, using a Docker image called nvidia-smi, which is a NVidia utility allowing to monitor (and manage) GPUs:

docker run --runtime=nvidia --rm nvidia/cuda nvidia-smi

Launching the previous command should return the following output:

+---





