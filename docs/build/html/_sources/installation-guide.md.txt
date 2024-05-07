# The SPARCED Installation Guide for Absolute Beginners
_Written by Aurore Amrit_

Hi! 🌄

If you are new to SPARCED and wish to get a working environment setup on
Ubuntu, then you are at the right place!
This is a document I wrote as a summer intern at the Birtwistle lab to make the
process easier for newcomers like me 🙂

## Environment
I am running an **Ubuntu 22.04 LTS** virtual machine on **VirtualBox**.
This guide should work even if you are using another hypervisor than VirtualBox
or that you are running Ubuntu directly on your computer.
With a few arrangements, the described steps should also work for other
versions of Ubuntu or any Debian-based Linux distribution.

:warning: Make sure you have enough disk space (30 GB is a minimum) :warning:

Before starting, open a terminal and run the following commands to make sure
everything is up to date:
```bash
sudo apt-get update
sudo apt-get upgrade
```
:coin: **Tip:** Ubuntu's terminal has autocompletion, so if you don't know the
end of the name of a file (for example because of version numbers) while
entering a command in your terminal, just press ```[tab]``` and see if it
fills in correctly.

## VirtualBox Guest Additions
_If you are not using VirtualBox, skip this section._

In the menu at the top of the virtual machine you are running, go to
_Devices_>_Insert Guest Additions CD image_, then in your terminal run:
```bash
sudo apt install build-essential dkms linux-headers-generic # Necessary if you chose the minimal installation of Ubuntu
lsblk | grep "rom"
cd /media/{username}/VBox_GAs_{version number} # with {username} being your username and {version number} the version number
sudo ./VBoxLinuxAdditions.run
```
You will need to restart the VM after that. Don't forget to eject the CD! 😉

## OpenMPI
_If you are not going to use parallel computation on your computer,
skip this section._

If you want to run some parallel code on your computer (for example to debug on
your own machine some code intented to run parallely on Palmetto), then you
will need to install OpenMPI.

:coin: **Tip:** Starting with the OpenMPI installation will prevent you from
starting everything from scratch again in case of failure (as it happened to me).

First, download the latest stable version of OpenMPI from the
[openmpi.org](https://www.open-mpi.org//software/ompi/v4.1/) website. You want
the ```.tar.gz``` extension.
Then run the following commands in your terminal (some can take a few minutes
and be very verbose, so get a ☕):
```bash
mkdir openmpi
cd openmpi
cp ~/Downloads/openmpi-{version number}.tar.gz . # with {version number} being the version number
tar -xzvf openmpi-{version number}.tar.gz # with {version number} being the version number
cd openmpi-{version number} # with {version number} being the version number
./configure --prefix=$HOME/openmpi # do not add any flag related to C++ (cxx) as they are no longer supported
make install
export PATH=$HOME/openmpi/bin:$PATH
export LD_LIBRARY_PATH=$HOME/openmpi/lib:$LD_LIBRARY_PATH
```
You can check if the installation process worked using:
```bash
mpirun --version
```
If you get an error involving Fortran during this process and you find a way to
fix it without having to reinstall everything, please notify me! 🙏

## Git, GitHub & SSH
### Git
```bash
sudo apt install git-all
```
You can use the following commands to set your username and your user's email.
```bash
git config --global user.name {username}
git config --global user.email {email}
```
### SSH for GitHub
I assume that you already have an account on [GitHub](https://github.com/).
```bash
ssh-keygen -t ed25519 -C "{email}" # with {email} being your email for GitHub
eval "$(ssh-agent -s)"
ssh-add ~/.ssh/id_ed25519
```
:coin: **Tip:** Just press ```[enter]``` to use the default file in which to
save the key.

Print your SSH key in the terminal using:
```bash
cat ~/.ssh/id_ed25519.pub
```
Copy it. The in your GitHub settings, create a new SSH key and paste it.

Run the following command to test your SSH connexion:
```bash
ssh -T git@github.com
```

## Anaconda
Download the
[Anaconda installer for Linux](https://www.anaconda.com/products/distribution#linux),
then run the following commands in your terminal:
```bash
sudo apt-get install libgl1-mesa-glx libegl1-mesa libxrandr2 libxrandr2 libxss1 libxcursor1 libxcomposite1 libasound2 libxi6 libxtst6
bash ~/Downloads/Anaconda3-{version-number}-Linux-x86_64.sh # with {version-number} being your version number
```
Then follow the instructions of the installer. Make sure it initializes
Anaconda3 by running conda init (type "yes" when asked for).

Once the installation is over, you need to restart your terminal (close and
reopen it or type ```source ~/.bashrc```).

Conda will automatically activate your ```base``` environment when launching
the terminal. If you want to disable this behavior, enter:
```bash
conda config --set auto_activate_base False # set it according to your preferences
```
Finally, verify your installation using:
```bash
conda list
```
To make sure Anaconda is up to date, run:
```bash
conda update -all
```

### Environment
Create a new Anaconda environment using the following commands:
```bash
conda create -n sparced # Creates an environment named "sparced"
source activate sparced # Activates the "sparced" environment
```
Unless you decide to set it otherwise, you will have to manually activate the
"sparced" environment each time you reopen your terminal.

### Python Packages
```bash
conda install matplotlib pandas scipy
pip install python-libsbml==5.18.0
pip install -Iv antimony==2.12.0.1 # WARNING: antimony >= 2.13.0 doesn't work with SPARCED
```
### The Amici Package
```bash
sudo apt install libatlas-base-dev swig
pip install amici==0.11.12 # WARNING: newer versions don't work with SPARCED
```
You might get an error about the CBLAS library (this happens mostly on
Palmetto), to fix it run:
```bash
conda install -c conda-forge openblas
export BLAS_LIBS=-lopenblas
```

### The mpi4py Package
_If you are not going to use parallel computation, skip this section._
```bash
conda remove compilers # if the compilers package is missing then don't install it!
conda install -c forge mpi4py # if you encounter any dependency version failure, try downgrading to Python 3.11 by typing 'conda install python=3.11'
conda install -c conda-forge compilers
python -m pip install gmx_MMPBSA
```

## TODO: Docker
_If you are not going to run SPARCED inside the official Jupyter Notebook
container, skip this section._

## SPARCED 🎆
This is only a setup suggestion:
```bash
cd ~/Documents
mkdir birtwistle-lab ; cd birtwistle-lab
git clone --recursive ssh://git@github.com/birtwistlelab/SPARCED.git # The official SPARCED repository
cd ..
git clone --recursive ssh://git@github.com/{username}/SPARCED.git # with {username} being your username on GitHub, assuming that you already forked SPARCED
```

## Clean
Remove all unused packages that were installed by dependencies during the setup:
```bash
sudo apt-get autoremove
```

Congratulations! You now have a full setup of SPARCED! 🦠
