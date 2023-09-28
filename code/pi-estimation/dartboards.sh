
# If we pass the "clean" key, then we delete all the existing images.
if [ "$1" == "clean" ]; then
    rm -rf output/figures/dartboards/*.*g
    python dartboards.py;
fi

# Convert all the images to gifs.
convert -resize 35% -delay 8 -loop 0 -deconstruct output/figures/dartboards/*.jpg output/figures/dartboard.gif
