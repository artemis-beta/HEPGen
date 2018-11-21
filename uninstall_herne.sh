pip uninstall -y $@ Herne &&
git checkout master &&
pip install -y --upgrade $@ .
