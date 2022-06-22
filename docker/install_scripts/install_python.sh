#!/bin/bash
env_name="core";
apt-get update -y && DEBIAN_FRONTEND=noninteractive apt-get install -y \
    python3 \
    python3-pip \
    python3-dev \
    libpq-dev \
    libsasl2-dev \
    libldap2-dev \
    libssl-dev \
    iputils-ping \
    python3-venv