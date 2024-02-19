# CellML Model Verification

This project focuses on verifying CellML models that cannot be directly converted to Bond Graphs.

## Overview

The code is designed to read CellML files with specially annotated variables, encrypted to contain specific information for decoding in Python.

## Project Structure

- **Organization:** Methods and functions are structured across various Python files, each dedicated to specific tasks.
- **Main Branch:** The most current and reliable version of the code is always available in the main branch.

## Getting Started

1. **Dependencies:** Install required packages listed in `requirements.txt`.
2. **CellML file:** Change the path in "cellml_file_dir" and "cellml_file" in main.py to your desired CellML file path.
3. **Run the Code:** Execute `main.py` to initiate the verification process.

## How It Works

The verification task is orchestrated by the `main.py` file, which calls essential functions for the process.

## Note

Currently, annotations in CellML files are assumed to be variable IDs encrypted with specific information. Feel free to check the documentation for more details.
