# Adds Herne compatibility for easy configuration

git clone https://github.com/artemis-beta/Herne.git &&

pip install --upgrade $@ Herne/ &&

git checkout hepgen-herne &&

pip install --upgrade $@ .
