pip uninstall -y Herne &&
rm -rf Herne
git checkout master &&
pip install --upgrade $@ .
