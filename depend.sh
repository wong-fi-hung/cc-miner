#!/bin/bash
clear
apt update -y
apt upgrade -y
apt install git wget nano && apt-get install automake autoconf pkg-config libcurl4-openssl-dev libjansson-dev libssl-dev libgmp-dev zlib1g-dev make g++ cmake libuv1-dev libhwloc-dev && sudo apt-get install libcurl4-openssl-dev libssl-dev libjansson-dev automake autotools-dev build-essential && sudo apt-get install git build-essential cmake libuv1-dev libssl-dev libhwloc-dev && apt update && apt upgrade && apt-get install libpci-dev
