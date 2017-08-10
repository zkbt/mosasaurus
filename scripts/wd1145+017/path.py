# here's how to install stuff locally:
# pip install --install-option="/home/hdiamond/local/lib/python2.7/site-packages" package-name
# if it's one of those python setup.py install situations, add the tag --user
# run this script at the beginning of every python session


import sys
sys.path.append('/home/hdiamond/local/lib/python2.7/site-packages/')
sys.path.append('/h/mulan0/code/')
sys.path.append('/h/mulan0/code/mosasaurus')
