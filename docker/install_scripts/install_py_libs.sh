#!/bin/bash
python3 -m pip install \
	redis \
	pika \
	numpy \
	fastapi \
	uvicorn[standard] \
	laspy[lazrs] \
	protobuf==3.20.* \
	h5py \
	python-dotenv \
	argparse