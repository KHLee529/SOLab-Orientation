# 10 Bar Truss Optimization

This repo contains codes for basic FEM (finite element method) and the using example of 10 bar truss optimization as the homework of SOLab orientation.

## Installation & Execution

Run the following commands to clone this repo to local storage.

```shell
$ git clone https://github.com/KHLee529/SOLab-Orientation.git
$ cd SOLab-Orientation
```

### Dependencies installation

Pipenv is used as the virtual environment and package management for this repo. \
It is recommended to use for isolating develop environment from other project.

Still, `requirements.txt` is provided for ones using pip or other virtual environment for package management.

- **pipenv**
	```shell
	$ pipenv install
	```
- **pip**
	```shell
	$ pip install -r requirements.txt
	```

### Execution

Run the following commands to do the 10 bar truss optimization. The optimization result will be printed as the final result.

```shell
$ python main.py
```

