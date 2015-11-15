GPM DNA Program README

Contents:
1. Introduction
2. Hardware and Software Requirements
3. Files
4. Binary install
5. Build from source
6. Contacts

1. Introduction

GPMDNA is a forensic DNA profile comparison program. The characteristic feature of this software is that genotype values at each loci are represented by genotype probability matrices (GPMs) allowing the storing and searching of electropherogram (epg) interpretations from both high and poor quality samples.

This distribution can be obtained from a statistical model, or from  trained human experts, thus extracting much information that would be discarded by standard interpretation rules. 

The program compares profiles arising from two DNA samples, one or both of which can have multiple donors and be affected by low DNA template or degraded DNA. It computes likelihood ratios to evaluate the hypothesis that the two samples have a common DNA donor, or two donors with a specified relatedness. It includes the standard population models, user-set FST values and a mutation model.

The program, particularly with optimised hardware, is capable of rapid searching of single or multiple profiles against large databases. This speed allows analysts an interactive capability particularly useful during complex, multi scene investigations.  The program can assist both in rapid identification of intelligence leads and in evaluating evidence within a judicial process. New and historical epgs that could not be used conventionally can now be searched and stored to the advantage of ongoing and "cold case" investigations.  We envisage crime datasets made up of probabilistic epg encodings that would be a subset of, or complementary to, conventional national crime databases.

GPMDNA is released under Crown copyright and the MIT open source licence to encourage discussion and acceptance of the probabalistic approach to forensic DNA matching in particular and to evidence in general. It is hoped that acceptance will lead to greater use of the approach in legal settings.

The software and supporting papers should reduce barriers to participation and support the idea of open standards and mutually beneficial development.

The method behind the GPMDNA program is described in detail in the following preprint:  XXX

2. Hardware and software requirements

The software compiles and runs on CentOS 5 linux. It should not be too hard to get it to work on other linux systems but this has not been tested. There are currently no windows/mac ports.

The software should run on any 64-bit processor. CUDA acceleration with NVIDIA GPUs is available, but this is optional. If there is no CUDA hardware then all calculations will automatically be run on the CPU. (NB you still need CUDA libraries, so the software can query them to find out you have no CUDA hardware! The binary install includes CUDA libraries). CUDA acceleration is advantageous only when matching two large sets of profiles against each other (or one large set of profiles against itself). For evaluation purposes, it is recommended that CPU matching is used.

3. Files

INSTALL.pdf:       Installation notes
LICENSE.txt:       License
README:            This README
binary:            Contains a binary install for CentOS 5.
build:             Freely available libraries and drivers needed to build the software
docs:              User manual and install notes
eclipse_workspace: Source code and project files for Eclipse Galileo
scripts:           Update script for building a binary install

4. Binary install

The binary install includes all required libraries, but you will need to create a database in MySQL and create/copy an allele frequency database. See the end of the INSTALL.pdf file for details. Edit the gmatch.sh script to suit your system.

5. Build from source

The way we build the software is with Eclipse Galileo. If it does not work
first time, try it again.

6. Contacts. 

For questions relating to copyright and licensing: keith2012roving@gmail.com

For questions relating to the software: gareth@dgwsoft.co.uk 

