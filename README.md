# mosasaurus
Tools for extracting chromatic lightcurves from MultiObject Spectrograph data. This

*These tools are still under active development. There's no promise they'll do something useful for you if you download them.*

### Installation


We recommend you fork/clone this repository onto your own computer and install directly from that editable package. From some base directory where you like to store code, run
```
git clone https://github.com/zkbt/mosasaurus.git
cd mosasaurus/
pip install -e .
```
This will download the current `mosasaurus` repository to your current directory, and then link the installed version of the `mosasaurus` package to this local repository. Changes you make to the code in the repository should be reflected in the version Python sees when it tries to `import mosasaurus`.

`mosasaurus` depends on [`craftroom`](https://github.com/zkbt/craftroom), a package of a few random astronomy tools. It's not on PyPI yet, so please run

`pip install git+https://github.com/zkbt/craftroom.git`

to install the latest. If you want to poke around inside `craftroom` too, you can use the same procedure to pip install it as your own editable local package repository as above.


### Usage

*(coming soon?)*

### Contributers

+ [Zach Berta-Thompson](http://casa.colorado.edu/~bertathompson/)
+ [Hannah Diamond-Lowe](https://www.cfa.harvard.edu/~hdiamondlowe/)
