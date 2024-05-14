#!/bin/bash

phonopy -f disp-001/vasprun.xml
phonopy --writefc --dim="6 6 6" --full-fc
