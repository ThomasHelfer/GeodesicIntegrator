for f in Image*.ppm; do
  convert ./"$f" ./"${f%.ppm}.png"
done
