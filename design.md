---
layout: default 
title: Design Guide
nav_order: 4
---

# SHEPHERD Design Guide

## Language Choice
This is an area of active debate. Proposed languages are:
* Julia
* Fortran
* C
* C++
* Rust

Each has their own strengths and weaknesses that would be highly beneficial for ths project.  

Selection criteria include:
* Ability to compile and deploy on remote (HPC) compute resources (all except Julia)
* Ability to write the fastest code (all)
* Ability to script and run models interactively  (Julia)
* Interoperability with existing sparse linear algebra libraries (C, C++, Fortran)
  * Intel MKL
  * AOCL Sparse
* Wide user base and developer skill (all except Rust?)