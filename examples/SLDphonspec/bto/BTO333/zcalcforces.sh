#!/bin/bash

phonopy -f disp-{001..003}/vasprun.xml
phonopy --writefc --dim="3 3 3" --full-fc
