# Multi-echo SMS-EPI sequence

## Usage

To create the sequence files, do:

1. Install Python environment (see below)
2. In MATLAB:
   ```matlab
   >> setup;
   >> main;

   ```


## Setup 

### Python setup
The CAIPIA ky-kz sampling pattern is generated using Rüdiger Stirnberg's Python code.

1. Get the code: 
Download from 
https://github.com/HarmonizedMRI/3DEPI.
For example, on the Linux command line, do:
    ```
    $ git clone git@github.com:HarmonizedMRI/3DEPI.git
    ```

2. Set up python virtual environment (recommended) and install required libraries;
    ```
    sudo apt install python3.13-venv
    python3 -m venv myvenv
    source myvenv/bin/activate
    pip install matplotlib
    pip install scipy
    ```

### MATLAB setup

Run setup.m


