pip uninstall -y $@ Herne &&
git checkout master &&
pip install --upgrade $@ .
