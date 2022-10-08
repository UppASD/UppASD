#!/bin/bash
ffmpeg -r 25 -i snap%05d.png -c:v libx265 -crf 28 -pix_fmt yuv420p -tag:v hvc1  output.mp4
