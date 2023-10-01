
# If we pass the "fresh" key, then we delete all the existing images --- a fresh
# start.
if [ "$1" == "fresh" ]; then
    rm -rf output/figures/dartboards/*.*g
    python dartboards.py;
fi

# Convert all the images to gifs.
cd output/figures/
# convert -delay 10 -loop 0 -deconstruct dartboards/*.jpg _dartboards.gif
gifsicle --resize 800x_ --colors 16 _dartboards.gif > dartboards.gif

