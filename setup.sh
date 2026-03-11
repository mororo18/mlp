#!/bin/bash

set -e

echo "=== MLP Development Environment Setup ==="

# Update package lists
echo "Updating package lists..."
sudo apt update

# Install all required tools
echo "Installing compilers and interpreters..."
sudo apt install -y \
    gcc \
    g++ \
    gfortran \
    make \
    golang \
    julia \
    lua5.3 \
    luajit \
    octave \
    pypy3 \
    valgrind \
    openjdk-21-jdk \
    nodejs

# .NET SDK (for csharp)
sudo apt install -y dotnet-sdk-8.0

# Python virtual environment
echo "Setting up Python virtual environment..."
python3 -m venv .venv
.venv/bin/pip install pandas psutil

echo "=== Setup complete! ==="
echo ""
echo "To activate Python environment: source .venv/bin/activate"
